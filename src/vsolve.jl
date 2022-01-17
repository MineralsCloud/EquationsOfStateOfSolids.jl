using PolynomialRoots: roots
using Roots: Order2, Newton, newton, find_zeros, find_zero

using .FiniteStrains: FromEulerianStrain, FromNaturalStrain

export vsolve

function vsolve(
    eos::PressureEquation{<:Murnaghan1st},
    p;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
)
    @unpack v0, b0, b′0 = getparam(eos)
    soln = [v0 * (1 + b′0 / b0 * p)^(-1 / b′0)]
    return _clamp(soln, bounds)
end
function vsolve(
    eos::PressureEquation{<:Murnaghan2nd},
    p;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    h = sqrt(2b0 * b″0 - b′0^2)
    k = b″0 * p + b′0
    numerator = exp(-2 / h * atan(p * h / (2b0 + p * b′0))) * v0
    denominator = (abs((k - h) / (k + h) * (b′0 + h) / (b′0 - h)))^(1 / h)
    soln = [numerator / denominator]
    return _clamp(soln, bounds)
end
function vsolve(
    eos::EnergyEquation{<:BirchMurnaghan2nd},
    e;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
)
    @unpack v0, b0, e0 = getparam(eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        soln = map(FromEulerianStrain(v0), [f, -f])
        return _clamp(soln, bounds)
    elseif Δ < 0
        return typeof(v0)[]  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function vsolve(
    eos::PressureEquation{<:BirchMurnaghan2nd},
    p;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos)
    # Solve f for (3 B0 f (2f + 1))^2 == p^2
    fs = roots(
        [-(p / 3b0)^2, 0, 1, 10, 40, 80, 80, 32];
        polish = true,
        epsilon = stopping_criterion,
    )
    soln = _strain2volume(eos, v0, fs, p, chop, rtol)
    return _clamp(soln, bounds)
end
function vsolve(
    eos::BulkModulusEquation{<:BirchMurnaghan2nd},
    b;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos)
    # Solve f for ((7f + 1) * (2f + 1)^(5/2))^2 == (b/b0)^2
    fs = roots(
        [1 - (b / b0)^2, 24, 229, 1130, 3160, 5072, 4368, 1568];
        polish = true,
        epsilon = stopping_criterion,
    )
    soln = _strain2volume(eos, v0, fs, b, chop, rtol)
    return _clamp(soln, bounds)
end
function vsolve(
    eos::EnergyEquation{<:BirchMurnaghan3rd},
    e;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    if b′0 == 4
        @warn "`b′0 == 4` for a `BirchMurnaghan3rd` is just a `BirchMurnaghan2nd`!"
        return vsolve(EnergyEquation(BirchMurnaghan2nd(v0, b0, e0)), e; bounds = bounds)
    else
        a = b′0 - 4
        r = 1 / 3a  # b = 1
        d = (e0 - e) / (9 / 2 * b0 * v0)
        p, q = -r^3 - d / 2a, -r^2
        Δ = p^2 + q^3
        fs = -r .+ if Δ > 0
            [cbrt(p + √Δ) + cbrt(p - √Δ)]  # Only 1 real solution
        elseif Δ < 0
            SIN, COS = sincos(acos(p / abs(r)^3) / 3)
            [2COS, -COS - √3 * SIN, -COS + √3 * SIN] .* abs(r)  # Verified
        elseif Δ == 0
            if p == q == 0
                [0]  # 3 reals are equal
            else  # p == -q != 0
                [2cbrt(p), -cbrt(p)]  # 2 real roots are equal, leaving 2 solutions
            end
        else
            @assert false "Δ == p^2 + q^3 == $Δ. this should never happen!"
        end  # solutions are strains
        soln = _strain2volume(eos, v0, fs, e, chop, rtol)
        return _clamp(soln, bounds)
    end
end
function vsolve(
    eos::PressureEquation{<:BirchMurnaghan3rd},
    p;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0, b′0 = getparam(eos)
    # Solve f for (f (2f + 1)^(5/2) [2 + 3f (b′0 - 4)])^2 - (p / (3b0/2))^2 = 0
    fs = roots(
        [
            -(2p / 3 / b0)^2,
            0,
            4,
            12b′0 - 8,
            9b′0^2 + 48b′0 - 176,
            90b′0^2 - 80(2 + 3b′0),
            360b′0^2 + 320(7 - 6b′0),
            720b′0^2 + 64(122 - 75b′0),
            720b′0^2 + 768(13 - 7b′0),
            288b′0^2 + 2304(2 - b′0),
        ];
        polish = true,
        epsilon = stopping_criterion,
    )
    soln = _strain2volume(eos, v0, fs, p, chop, rtol)
    return _clamp(soln, bounds)
end
function vsolve(
    eos::EnergyEquation{<:BirchMurnaghan4th},
    e;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    h = b0 * b″0 + b′0^2
    fs = roots(
        [(e0 - e) / (3 / 8 * v0 * b0), 0, 12, 12(b′0 - 4), 143 - 63b′0 + 9h];
        polish = true,
        epsilon = stopping_criterion,
    )
    soln = _strain2volume(eos, v0, fs, e, chop, rtol)
    return _clamp(soln, bounds)
end
function vsolve(
    eos::EnergyEquation{<:PoirierTarantola2nd},
    e;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
)
    @unpack v0, b0, e0 = getparam(eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        soln = map(FromNaturalStrain(v0), [f, -f])
        return _clamp(soln, bounds)
    elseif Δ < 0
        return typeof(v0)[]  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function vsolve(
    eos::EquationOfStateOfSolids,
    y;
    bounds = (zero(eos.param.v0), 4 * eos.param.v0),
    xrtol = eps(),
    rtol = 4eps(),
)
    # Bisection method
    try
        return find_zeros(v -> eos(v) - y, bounds; xrtol = xrtol, rtol = rtol)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(eos.param.v0)[]
    end
end
function vsolve(
    eos::EquationOfStateOfSolids,
    y,
    vᵢ,
    method = Order2();
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    maxiter = 40,
    verbose = false,
    xrtol = eps(),
    rtol = 4eps(),
)
    try
        vᵣ = find_zero(
            v -> eos(v) - y,
            vᵢ,
            method;
            maxevals = maxiter,
            verbose = verbose,
            xrtol = xrtol,
            rtol = rtol,
        )
        return _clamp([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end
function vsolve(
    eos::EnergyEquation,
    e,
    vᵢ,
    ::Newton;
    bounds = (zero(eos.param.v0), Inf * eos.param.v0),
    kwargs...,
)
    try
        vᵣ = newton(v -> eos(v) - e, v -> -PressureEquation(eos)(v), vᵢ; kwargs...)
        return _clamp([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end

function _strain2volume(
    eos::EquationOfStateOfSolids{<:BirchMurnaghan},
    v0,
    fs,
    y,
    chop = eps(),
    rtol = sqrt(eps()),
)
    @assert _ispositive(chop) && _ispositive(rtol) "either `chop` or `rtol` is less than 0!"
    return fs |>
           Base.Fix1(map, FromEulerianStrain(v0)) |>
           Base.Fix1(filter, x -> abs(imag(x)) < chop * oneunit(imag(x))) .|>  # If `x` has unit
           real |>
           Base.Fix1(filter, x -> isapprox(eos(x), y; rtol = rtol)) |>  # In case of duplicate values
           unique
end

function _clamp(soln, bounds)
    low, high = extrema(bounds)
    return filter(v -> low <= v <= high, soln)
end
