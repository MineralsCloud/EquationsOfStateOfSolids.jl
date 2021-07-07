module Inverse

using Chain: @chain
using Configurations: from_kwargs, @option
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
    BirchMurnaghan4th,
    PoirierTarantola,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    getparam
using ..FiniteStrains: FromEulerianStrain, FromNaturalStrain

export NumericalInversionOptions

abstract type Inverted{T<:EquationOfStateOfSolids} end
struct AnalyticallyInverted{T} <: Inverted{T}
    eos::T
end
struct NumericallyInverted{T} <: Inverted{T}
    eos::T
end

@option "num_inv" struct NumericalInversionOptions
    search_interval::AbstractVector = [eps(), 2]
    maxiter::Int64 = 40
    verbose::Bool = false
end

function (x::AnalyticallyInverted{<:PressureEquation{<:Murnaghan1st}})(p)
    @unpack v0, b0, b′0 = getparam(x.eos)
    return [v0 * (1 + b′0 / b0 * p)^(-1 / b′0)]
end
function (x::AnalyticallyInverted{<:PressureEquation{<:Murnaghan2nd}})(p)
    @unpack v0, b0, b′0, b″0 = getparam(x.eos)
    h = sqrt(2b0 * b″0 - b′0^2)
    k = b″0 * p + b′0
    numerator = exp(-2 / h * atan(p * h / (2b0 + p * b′0))) * v0
    denominator = (abs((k - h) / (k + h) * (b′0 + h) / (b′0 - h)))^(1 / h)
    return [numerator / denominator]
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:BirchMurnaghan2nd}})(e)
    @unpack v0, b0, e0 = getparam(x.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromEulerianStrain(v0), [f, -f])
    elseif Δ < 0
        return []  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:BirchMurnaghan3rd}})(e)
    @unpack v0, b0, b′0, e0 = getparam(x.eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    if b′0 == 4
        @warn "`b′0 == 4` for a `BirchMurnaghan3rd` is just a `BirchMurnaghan2nd`!"
        eos⁻¹ = EnergyEquation(BirchMurnaghan2nd(v0, b0, e0))^(-1)
        return eos⁻¹(e)
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
        return @chain fs begin
            map(FromEulerianStrain(v0), _)
            filter(isreal, _)
            @. real
            filter(_ispositive, _)
        end
    end
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:BirchMurnaghan4th}})(e)
    @unpack v0, b0, b′0, b″0, e0 = getparam(x.eos)
    h = b0 * b″0 + b′0^2
    fs = roots([e0 - e, 3 // 8 * v0 * b0 .* (9h - 63b′0 + 143, 12 * (b′0 - 4), 12)...])
    return @chain fs begin
        map(FromEulerianStrain(v0), _)
        filter(isreal, _)
        @. real
        filter(_ispositive, _)
    end
end
function (x::AnalyticallyInverted{<:EnergyEquation{<:PoirierTarantola2nd}})(e)
    @unpack v0, b0, e0 = getparam(x.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromNaturalStrain(v0), [f, -f])
    elseif Δ < 0
        return []  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(
    y,
    method::Union{AbstractBracketing,AbstractSecant},
    options::NumericalInversionOptions,
)
    v0 = _within(options.search_interval, method) .* getparam(x.eos).v0  # v0 can be negative
    @assert _ispositive(minimum(v0))  # No negative volume
    v = find_zero(
        guess -> x.eos(guess) - y,
        v0,
        method;
        maxevals = options.maxiter,
        verbose = options.verbose,
    )
    if !_ispositive(v)
        @warn "the volume found is negative!"
    end
    return v
end
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(
    y,
    method::Union{AbstractBracketing,AbstractSecant};
    kwargs...,
)
    options = from_kwargs(NumericalInversionOptions; kwargs...)
    return x(y, method, options)
end
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(
    y,
    options::NumericalInversionOptions,
)
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
        if options.verbose
            @info "using method `$T`..."
        end
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return x(y, T(), options)
        catch e
            if options.verbose
                @info "method `$T` failed because of `$e`."
            end
            continue
        end
    end
    error("no volume found!")
end
function (x::NumericallyInverted{<:EquationOfStateOfSolids})(y; kwargs...)
    options = from_kwargs(NumericalInversionOptions; kwargs...)
    return x(y, options)
end
_within(search_interval, ::AbstractBracketing) = extrema(search_interval)
_within(search_interval, ::AbstractSecant) = sum(extrema(search_interval)) / 2

# Idea from https://discourse.julialang.org/t/functional-inverse/10959/6
Base.literal_pow(::typeof(^), eos::EquationOfStateOfSolids, ::Val{-1}) =
    NumericallyInverted(eos)
Base.literal_pow(
    ::typeof(^),
    eos::Union{
        PressureEquation{<:Murnaghan},
        EnergyEquation{<:BirchMurnaghan},
        EnergyEquation{<:PoirierTarantola},
    },
    ::Val{-1},
) = AnalyticallyInverted(eos)

end
