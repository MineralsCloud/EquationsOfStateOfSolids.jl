module Inverse

using Compat: filter
using InteractiveUtils: subtypes
using PolynomialRoots: roots
using Roots:
    find_zero,
    AbstractBracketing,
    AbstractSecant,
    Bisection,
    BisectionExact,
    FalsePosition,
    A42,
    AlefeldPotraShi,
    Brent,
    Esser,
    King,
    KumarSinghAkanksha,
    Order0,
    Order16,
    Order1B,
    Order2,
    Order2B,
    Order5,
    Order8,
    Secant,
    Steffensen,
    Thukral16,
    Thukral8
using UnPack: @unpack

using ..EquationsOfStateOfSolids:
    _ispositive,
    PressureEquation,
    EnergyEquation,
    EquationOfStateOfSolids,
    Murnaghan,
    Murnaghan1st,
    Murnaghan2nd,
    BirchMurnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    PoirierTarantola,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    getparam
using ..FiniteStrains: FromEulerianStrain, FromNaturalStrain

export inverse

abstract type Inverted{T<:EquationOfStateOfSolids} end
struct AnalyticallyInverted{T} <: Inverted{T}
    eos::T
end
struct NumericallyInverted{T} <: Inverted{T}
    eos::T
end

function (x::AnalyticallyInverted{<:PressureEquation{<:Murnaghan1st}})(p)
    @unpack v0, b0, b′0, e0 = getparam(x.eos)
    return (v0 * (1 + b′0 / b0 * p)^(-1 / b′0),)
end
function (x::AnalyticallyInverted{<:PressureEquation{<:Murnaghan2nd}})(p)
    @unpack v0, b0, b′0, b″0 = getparam(x.eos)
    h = sqrt(2b0 * b″0 - b′0^2)
    k = b″0 * p + b′0
    term1 = exp(-2 / h * atan(p * h / (2b0 + p * b′0))) * v0
    term2 = (abs((k - h) / (k + h) * (b′0 + h) / (b′0 - h)))^(1 / h)
    return term1 / term2
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:BirchMurnaghan2nd}})(e)
    @unpack v0, b0, b′0, e0 = getparam(x.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromEulerianStrain(v0), (f, -f))
    elseif Δ < 0
        return ()  # Complex strains
    else
        throw(ArgumentError("this should never happen!"))
    end
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:BirchMurnaghan3rd}})(e)
    @unpack v0, b0, b′0, e0 = getparam(x.eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    a = b′0 - 4
    r = 1 / 3a  # b = 1
    d = (e0 - e) / (9 / 2 * b0 * v0)
    p, q = -r^3 - d / 2a, -r^2
    Δ = p^2 + q^3
    fs = -r .+ if Δ > 0
        (cbrt(p + √Δ) + cbrt(p - √Δ),)  # Only 1 real solution
    elseif Δ < 0
        SIN, COS = sincos(acos(p / abs(r)^3) / 3)
        (2COS, -COS - √3 * SIN, -COS + √3 * SIN) .* abs(r)  # Verified
    elseif Δ == 0
        if p == q == 0
            (0,)  # 3 reals are equal
        else  # p == -q != 0
            2cbrt(p), -cbrt(p)  # 2 real roots are equal, leaving 2 solutions
        end
    else
        throw(ArgumentError("this should never happen!"))
    end  # solutions are strains
    vs = map(FromEulerianStrain(v0), fs)
    return filter(_ispositive, map(real, filter(isreal, vs)))
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:PoirierTarantola2nd}})(e)
    @unpack v0, b0, b′0, e0 = getparam(x.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromNaturalStrain(v0), (f, -f))
    elseif Δ < 0
        return ()  # Complex strains
    else
        throw(ArgumentError("this should never happen!"))
    end
end
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(
    y,
    method::Union{AbstractBracketing,AbstractSecant};
    interval = (0.5, 1.5),
    maxiter = 40,
    verbose = false,
)
    v0 = _within(interval, method) .* getparam(x.eos).v0  # v0 can be negative
    @assert _ispositive(minimum(v0))  # No negative volume
    v = find_zero(guess -> x.eos(guess) - y, v0, method; maxevals = maxiter, verbose = verbose)
    if !_ispositive(v)
        @warn "the volume found is negative!"
    end
    return v
end
_within(interval, ::AbstractBracketing) = extrema(interval)
_within(interval, ::AbstractSecant) = sum(extrema(interval)) / 2
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(y; verbose = false, kwargs...)
    for T in [
        Bisection,
        BisectionExact,
        FalsePosition,
        A42,
        AlefeldPotraShi,
        Brent,
        Esser,
        King,
        KumarSinghAkanksha,
        Order0,
        Order16,
        Order1B,
        Order2,
        Order2B,
        Order5,
        Order8,
        Secant,
        Steffensen,
        Thukral16,
        Thukral8,
    ]
        if verbose
            @info "using method `$T`..."
        end
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return x(y, T(); verbose = verbose, kwargs...)
        catch e
            if verbose
                @info "method `$T` failed because of `$e`."
            end
            continue
        end
    end
    error("no volume found!")
end

inverse(eos::EquationOfStateOfSolids) = NumericallyInverted(eos)
inverse(eos::PressureEquation{<:Murnaghan}) = AnalyticallyInverted(eos)
inverse(eos::EnergyEquation{<:BirchMurnaghan}) = AnalyticallyInverted(eos)
inverse(eos::EnergyEquation{<:PoirierTarantola}) = AnalyticallyInverted(eos)

end
