using PolynomialRoots: roots
using Roots: Order2, Newton, newton, find_zeros, find_zero

using .FiniteStrains: FromEulerianStrain, FromNaturalStrain

abstract type VolumeSolver end
struct AnalyticalVolumeSolver{E<:EquationOfStateOfSolids,V} <: VolumeSolver
    eos::E
    bounds::NTuple{2,V}
end
function AnalyticalVolumeSolver(eos, bounds=(0, Inf) .* eos.param.v0)
    return AnalyticalVolumeSolver(eos, extrema(bounds))
end
function EnergySolver(params::Parameters, bounds=(0, Inf) .* params.v0)
    return AnalyticalVolumeSolver(EnergyEquation(params), bounds)
end
function PressureSolver(params::Parameters, bounds=(0, Inf) .* params.v0)
    return AnalyticalVolumeSolver(PressureEquation(params), bounds)
end
function BulkModulusSolver(params::Parameters, bounds=(0, Inf) .* params.v0)
    return AnalyticalVolumeSolver(BulkModulusEquation(params), bounds)
end

struct NumericalVolumeSolver{E<:EquationOfStateOfSolids,M,V} <: VolumeSolver
    eos::E
    method::M
    bounds::NTuple{2,V}
end
function NumericalVolumeSolver(eos, method, bounds=(0, Inf) .* eos.param.v0)
    return NumericalVolumeSolver(eos, method, extrema(bounds))
end
function EnergySolver(params::Parameters, method, bounds=(0, Inf) .* params.v0)
    return NumericalVolumeSolver(EnergyEquation(params), method, bounds)
end
function PressureSolver(params::Parameters, method, bounds=(0, Inf) .* params.v0)
    return NumericalVolumeSolver(PressureEquation(params), method, bounds)
end
function BulkModulusSolver(params::Parameters, method, bounds=(0, Inf) .* params.v0)
    return NumericalVolumeSolver(BulkModulusEquation(params), method, bounds)
end

function (::AnalyticalVolumeSolver{<:PressureEquation{<:Murnaghan1st}})(p)
    eos, bounds = problem.eos, problem.bounds
    @unpack v0, b0, b′0 = getparam(eos)
    solution = [v0 * (1 + b′0 / b0 * p)^(-1 / b′0)]
    return sieve(solution, bounds)
end
function (::AnalyticalVolumeSolver{<:PressureEquation{<:Murnaghan2nd}})(p)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    h = sqrt(2b0 * b″0 - b′0^2)
    k = b″0 * p + b′0
    numerator = exp(-2 / h * atan(p * h / (2b0 + p * b′0))) * v0
    denominator = (abs((k - h) / (k + h) * (b′0 + h) / (b′0 - h)))^(1 / h)
    solution = [numerator / denominator]
    return sieve(solution, bounds)
end
function (::AnalyticalVolumeSolver{<:EnergyEquation{<:BirchMurnaghan2nd}})(e)
    @unpack v0, b0, e0 = getparam(eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        solutions = map(FromEulerianStrain(v0), [f, -f])
        return sieve(solutions, bounds)
    elseif Δ < 0
        return typeof(v0)[]  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (::AnalyticalVolumeSolver{<:PressureEquation{<:BirchMurnaghan2nd}})(
    p;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos)
    # Solve f for (3 B0 f (2f + 1))^2 == p^2
    fs = roots(
        [-(p / 3b0)^2, 0, 1, 10, 40, 80, 80, 32]; polish=true, epsilon=stopping_criterion
    )
    solutions = _strain2volume(eos, v0, fs, p, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalVolumeSolver{<:BulkModulusEquation{<:BirchMurnaghan2nd}},
    b;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos)
    # Solve f for ((7f + 1) * (2f + 1)^(5/2))^2 == (b/b0)^2
    fs = roots(
        [1 - (b / b0)^2, 24, 229, 1130, 3160, 5072, 4368, 1568];
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v0, fs, b, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalVolumeSolver{<:EnergyEquation{<:BirchMurnaghan3rd}},
    e;
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    if b′0 == 4
        @warn "`b′0 == 4` for a `BirchMurnaghan3rd` is just a `BirchMurnaghan2nd`!"
        return solve(EnergyEquation(BirchMurnaghan2nd(v0, b0, e0)), e; bounds=bounds)
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
        solutions = _strain2volume(eos, v0, fs, e, chop, rtol)
        return sieve(solutions, bounds)
    end
end
function (
    ::AnalyticalVolumeSolver{<:PressureEquation{<:BirchMurnaghan3rd}},
    p;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
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
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v0, fs, p, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalVolumeSolver{<:EnergyEquation{<:BirchMurnaghan4th}},
    e;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    h = b0 * b″0 + b′0^2
    fs = roots(
        [(e0 - e) / (3 / 8 * v0 * b0), 0, 12, 12(b′0 - 4), 143 - 63b′0 + 9h];
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v0, fs, e, chop, rtol)
    return sieve(solutions, bounds)
end
function (::AnalyticalVolumeSolver{<:EnergyEquation{<:PoirierTarantola2nd}})(e)
    @unpack v0, b0, e0 = getparam(eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        solutions = map(FromNaturalStrain(v0), [f, -f])
        return sieve(solutions, bounds)
    elseif Δ < 0
        return typeof(v0)[]  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (::AnalyticalVolumeSolver)(y; xrtol=eps(), rtol=4eps())
    # Bisection method
    try
        return find_zeros(v -> eos(v) - y, bounds; xrtol=xrtol, rtol=rtol)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(eos.param.v0)[]
    end
end
function (::NumericalVolumeSolver{E,Order2})(
    y, vᵢ; maxiter=40, verbose=false, xrtol=eps(), rtol=4eps()
) where {E}
    try
        vᵣ = find_zero(
            v -> eos(v) - y,
            vᵢ,
            method;
            maxevals=maxiter,
            verbose=verbose,
            xrtol=xrtol,
            rtol=rtol,
        )
        return sieve([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end
function (::NumericalVolumeSolver{E,Newton})(e, vᵢ; kwargs...) where {E}
    try
        vᵣ = newton(v -> eos(v) - e, v -> -PressureEquation(eos)(v), vᵢ; kwargs...)
        return sieve([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end

function _strain2volume(
    eos::EquationOfStateOfSolids{<:BirchMurnaghan}, v0, fs, y, chop=eps(), rtol=sqrt(eps())
)
    @assert !isnegative(chop) && !isnegative(rtol) "either `chop` or `rtol` is less than 0!"
    return unique(
        Base.Fix1(filter, x -> isapprox(eos(x), y; rtol=rtol))(
            real.(
                Base.Fix1(filter, x -> abs(imag(x)) < chop * oneunit(imag(x)))(
                    Base.Fix1(map, FromEulerianStrain(v0))(fs)
                )
            ),
        ),
    )
end

function sieve(solutions, bounds)
    lower, upper = extrema(bounds)
    return filter(v -> lower <= v <= upper, solutions)
end
