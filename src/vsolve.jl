using PolynomialRoots: roots
using Roots:
    Order2,
    AbstractNewtonLikeMethod,
    AbstractHalleyLikeMethod,
    AbstractSecantMethod,
    LithBoonkkampIJzerman,
    find_zeros,
    find_zero

using .FiniteStrains: VolumeFromEulerianStrain, VolumeFromNaturalStrain

abstract type VolumeSolver end
struct AnalyticalSolver{E<:EquationOfStateOfSolids,V} <: VolumeSolver
    eos::E
    bounds::NTuple{2,V}
end
function AnalyticalSolver(eos, bounds=(0, Inf) .* eos.param.v₀)
    return AnalyticalSolver(eos, extrema(bounds))
end
function EnergySolver(params::Parameters, bounds=(0, Inf) .* params.v₀)
    return AnalyticalSolver(EnergyEquation(params), bounds)
end
function PressureSolver(params::Parameters, bounds=(0, Inf) .* params.v₀)
    return AnalyticalSolver(PressureEquation(params), bounds)
end
function BulkModulusSolver(params::Parameters, bounds=(0, Inf) .* params.v₀)
    return AnalyticalSolver(BulkModulusEquation(params), bounds)
end

struct NumericalSolver{E<:EquationOfStateOfSolids,M,V} <: VolumeSolver
    eos::E
    method::M
    bounds::NTuple{2,V}
end
function NumericalSolver(eos, method, bounds=(0, Inf) .* eos.param.v₀)
    return NumericalSolver(eos, method, extrema(bounds))
end
function EnergySolver(params::Parameters, method, bounds=(0, Inf) .* params.v₀)
    return NumericalSolver(EnergyEquation(params), method, bounds)
end
function PressureSolver(params::Parameters, method, bounds=(0, Inf) .* params.v₀)
    return NumericalSolver(PressureEquation(params), method, bounds)
end
function BulkModulusSolver(params::Parameters, method, bounds=(0, Inf) .* params.v₀)
    return NumericalSolver(BulkModulusEquation(params), method, bounds)
end

function (::AnalyticalSolver{<:PressureEquation{<:Murnaghan1st}})(p)
    eos, bounds = problem.eos, problem.bounds
    v₀, b₀, b′₀ = unpack(eos)
    solution = [v₀ * (1 + b′₀ / b₀ * p)^(-1 / b′₀)]
    return sieve(solution, bounds)
end
function (::AnalyticalSolver{<:PressureEquation{<:Murnaghan2nd}})(p)
    v₀, b₀, b′₀, b″₀ = unpack(eos)
    h = sqrt(2b₀ * b″₀ - b′₀^2)
    k = b″₀ * p + b′₀
    numerator = exp(-2 / h * atan(p * h / (2b₀ + p * b′₀))) * v₀
    denominator = (abs((k - h) / (k + h) * (b′₀ + h) / (b′₀ - h)))^(1 / h)
    solution = [numerator / denominator]
    return sieve(solution, bounds)
end
function (::AnalyticalSolver{<:EnergyEquation{<:BirchMurnaghan2nd}})(e)
    v₀, b₀, e₀ = unpack(eos)
    Δ = (e - e₀) / v₀ / b₀
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        solutions = map(VolumeFromEulerianStrain(v₀), [f, -f])
        return sieve(solutions, bounds)
    elseif Δ < 0
        return typeof(v₀)[]  # Complex strains
    else
        @assert false "Δ == (e - e₀) / v₀ / b₀ == $Δ. this should never happen!"
    end
end
function (::AnalyticalSolver{<:PressureEquation{<:BirchMurnaghan2nd}})(
    p;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v₀, b₀ = eos
    # Solve f for (3 B0 f (2f + 1))^2 == p^2
    fs = roots(
        [-(p / 3b₀)^2, 0, 1, 10, 40, 80, 80, 32]; polish=true, epsilon=stopping_criterion
    )
    solutions = _strain2volume(eos, v₀, fs, p, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalSolver{<:BulkModulusEquation{<:BirchMurnaghan2nd}},
    b;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v₀, b₀ = eos
    # Solve f for ((7f + 1) * (2f + 1)^(5/2))^2 == (b/b₀)^2
    fs = roots(
        [1 - (b / b₀)^2, 24, 229, 1130, 3160, 5072, 4368, 1568];
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v₀, fs, b, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalSolver{<:EnergyEquation{<:BirchMurnaghan3rd}},
    e;
    chop=eps(),
    rtol=sqrt(eps()),
)
    v₀, b₀, b′₀, e₀ = unpack(eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    if b′₀ == 4
        @warn "`b′₀ == 4` for a `BirchMurnaghan3rd` is just a `BirchMurnaghan2nd`!"
        return solve(EnergyEquation(BirchMurnaghan2nd(v₀, b₀, e₀)), e; bounds=bounds)
    else
        a = b′₀ - 4
        r = 1 / 3a  # b = 1
        d = (e₀ - e) / (9 / 2 * b₀ * v₀)
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
        solutions = _strain2volume(eos, v₀, fs, e, chop, rtol)
        return sieve(solutions, bounds)
    end
end
function (
    ::AnalyticalSolver{<:PressureEquation{<:BirchMurnaghan3rd}},
    p;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    v₀, b₀, b′₀ = unpack(eos)
    # Solve f for (f (2f + 1)^(5/2) [2 + 3f (b′₀ - 4)])^2 - (p / (3b₀/2))^2 = 0
    fs = roots(
        [
            -(2p / 3 / b₀)^2,
            0,
            4,
            12b′₀ - 8,
            9b′₀^2 + 48b′₀ - 176,
            90b′₀^2 - 80(2 + 3b′₀),
            360b′₀^2 + 320(7 - 6b′₀),
            720b′₀^2 + 64(122 - 75b′₀),
            720b′₀^2 + 768(13 - 7b′₀),
            288b′₀^2 + 2304(2 - b′₀),
        ];
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v₀, fs, p, chop, rtol)
    return sieve(solutions, bounds)
end
function (
    ::AnalyticalSolver{<:EnergyEquation{<:BirchMurnaghan4th}},
    e;
    stopping_criterion=1e-20,  # Unitless
    chop=eps(),
    rtol=sqrt(eps()),
)
    @unpack v₀, b₀, b′₀, b″₀, e₀ = eos
    h = b₀ * b″₀ + b′₀^2
    fs = roots(
        [(e₀ - e) / (3 / 8 * v₀ * b₀), 0, 12, 12(b′₀ - 4), 143 - 63b′₀ + 9h];
        polish=true,
        epsilon=stopping_criterion,
    )
    solutions = _strain2volume(eos, v₀, fs, e, chop, rtol)
    return sieve(solutions, bounds)
end
function (::AnalyticalSolver{<:EnergyEquation{<:PoirierTarantola2nd}})(e)
    v₀, b₀, e₀ = unpack(eos)
    Δ = (e - e₀) / v₀ / b₀
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        solutions = map(VolumeFromNaturalStrain(v₀), [f, -f])
        return sieve(solutions, bounds)
    elseif Δ < 0
        return typeof(v₀)[]  # Complex strains
    else
        @assert false "Δ == (e - e₀) / v₀ / b₀ == $Δ. this should never happen!"
    end
end
function (::NumericalSolver)(y; xrtol=eps(), rtol=4eps())
    # Bisection method
    try
        return find_zeros(v -> eos(v) - y, bounds; xrtol=xrtol, rtol=rtol)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(eos.param.v₀)[]
    end
end
function (
    ::NumericalSolver{
        <:EnergyEquation,<:Union{AbstractNewtonLikeMethod,LithBoonkkampIJzerman{1,1}}
    }
)(
    e, vᵢ; kwargs...
)  # Newton's method
    try
        vᵣ = find_zero(
            (v -> eos(v) - e, v -> -PressureEquation(eos)(v)), vᵢ, method; kwargs...
        )
        return sieve([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end
function (
    ::NumericalSolver{
        <:EnergyEquation,<:Union{AbstractHalleyLikeMethod,LithBoonkkampIJzerman{1,2}}
    }
)(
    e, vᵢ; kwargs...
)  # Halley-like method
    try
        vᵣ = find_zero(
            (
                v -> eos(v) - e,  # f
                v -> -PressureEquation(eos)(v),  # f′
                v -> BulkModulusEquation(eos)(v) / v,  # f″
            ),
            vᵢ,
            method;
            kwargs...,
        )
        return sieve([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end
function (
    ::NumericalSolver{
        <:EnergyEquation,<:Union{AbstractSecantMethod,LithBoonkkampIJzerman{2,0}}
    }
)(
    e, vᵢ; kwargs...
)  # The secant method
    try
        vᵣ = find_zero(v -> eos(v) - e, vᵢ, method; kwargs...)
        return sieve([vᵣ], bounds)
    catch e
        @error "cannot find solution! Come across `$e`!"
        return typeof(vᵢ)[]
    end
end

function _strain2volume(
    eos::EquationOfStateOfSolids{<:BirchMurnaghan}, v₀, fs, y, chop=eps(), rtol=sqrt(eps())
)
    @assert !isnegative(chop) && !isnegative(rtol) "either `chop` or `rtol` is less than 0!"
    return unique(
        Base.Fix1(filter, x -> isapprox(eos(x), y; rtol=rtol))(
            real.(
                Base.Fix1(filter, x -> abs(imag(x)) < chop * oneunit(imag(x)))(
                    Base.Fix1(map, VolumeFromEulerianStrain(v₀))(fs)
                )
            ),
        ),
    )
end

function sieve(solutions, bounds)
    lower, upper = extrema(bounds)
    return filter(v -> lower <= v <= upper, solutions)
end
