module Inverse

using PolynomialRoots: roots
using Roots: find_zero, Order16
using StaticArrays: MVector
using UnPack: @unpack

using ..EquationsOfStateOfSolids:
    _ispositive,
    BulkModulusEquation,
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

"Wrap an equation of state for inversions."
struct Inverted{T<:EquationOfStateOfSolids}
    eos::T
end

function (eos⁻¹::Inverted{<:PressureEquation{<:Murnaghan1st}})(p)
    @unpack v0, b0, b′0 = getparam(eos⁻¹.eos)
    return MVector(v0 * (1 + b′0 / b0 * p)^(-1 / b′0))
end
function (eos⁻¹::Inverted{<:PressureEquation{<:Murnaghan2nd}})(p)
    @unpack v0, b0, b′0, b″0 = getparam(eos⁻¹.eos)
    h = sqrt(2b0 * b″0 - b′0^2)
    k = b″0 * p + b′0
    numerator = exp(-2 / h * atan(p * h / (2b0 + p * b′0))) * v0
    denominator = (abs((k - h) / (k + h) * (b′0 + h) / (b′0 - h)))^(1 / h)
    return MVector(numerator / denominator)
end
function (eos⁻¹::Inverted{<:EnergyEquation{<:BirchMurnaghan2nd}})(e)
    @unpack v0, b0, e0 = getparam(eos⁻¹.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromEulerianStrain(v0), MVector(f, -f))
    elseif Δ < 0
        return MVector{0,typeof(v0)}()  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (eos⁻¹::Inverted{<:PressureEquation{<:BirchMurnaghan2nd}})(
    p;
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos⁻¹.eos)
    # Solve f for (3 B0 f (2f + 1))^2 == p^2
    fs = roots(
        [-(p / 3b0)^2, 0, 1, 10, 40, 80, 80, 32];
        polish = true,
        epsilon = stopping_criterion,
    )
    return _strain2volume(eos⁻¹.eos, v0, fs, p, chop, rtol)
end
function (eos⁻¹::Inverted{<:BulkModulusEquation{<:BirchMurnaghan2nd}})(
    b;
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0 = getparam(eos⁻¹.eos)
    # Solve f for ((7f + 1) * (2f + 1)^(5/2))^2 == (b/b0)^2
    fs = roots(
        [1 - (b / b0)^2, 24, 229, 1130, 3160, 5072, 4368, 1568];
        polish = true,
        epsilon = stopping_criterion,
    )
    return _strain2volume(eos⁻¹.eos, v0, fs, b, chop, rtol)
end
function (eos⁻¹::Inverted{<:EnergyEquation{<:BirchMurnaghan3rd}})(
    e;
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0, b′0, e0 = getparam(eos⁻¹.eos)
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
            MVector(cbrt(p + √Δ) + cbrt(p - √Δ))  # Only 1 real solution
        elseif Δ < 0
            SIN, COS = sincos(acos(p / abs(r)^3) / 3)
            MVector(2COS, -COS - √3 * SIN, -COS + √3 * SIN) .* abs(r)  # Verified
        elseif Δ == 0
            if p == q == 0
                MVector(0)  # 3 reals are equal
            else  # p == -q != 0
                MVector(2cbrt(p), -cbrt(p))  # 2 real roots are equal, leaving 2 solutions
            end
        else
            @assert false "Δ == p^2 + q^3 == $Δ. this should never happen!"
        end  # solutions are strains
        return _strain2volume(eos⁻¹.eos, v0, fs, e, chop, rtol)
    end
end
function (eos⁻¹::Inverted{<:EnergyEquation{<:BirchMurnaghan4th}})(
    e;
    stopping_criterion = 1e-20,  # Unitless
    chop = eps(),
    rtol = sqrt(eps()),
)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos⁻¹.eos)
    h = b0 * b″0 + b′0^2
    fs = roots(
        [(e0 - e) / (3 / 8 * v0 * b0), 0, 12, 12(b′0 - 4), 143 - 63b′0 + 9h];
        polish = true,
        epsilon = stopping_criterion,
    )
    return _strain2volume(eos⁻¹.eos, v0, fs, e, chop, rtol)
end
function (eos⁻¹::Inverted{<:EnergyEquation{<:PoirierTarantola2nd}})(e)
    @unpack v0, b0, e0 = getparam(eos⁻¹.eos)
    Δ = (e - e0) / v0 / b0
    if Δ >= 0
        f = sqrt(2 / 9 * Δ)
        return map(FromNaturalStrain(v0), MVector(f, -f))
    elseif Δ < 0
        return MVector{0,typeof(v0)}()  # Complex strains
    else
        @assert false "Δ == (e - e0) / v0 / b0 == $Δ. this should never happen!"
    end
end
function (eos⁻¹::Inverted{<:EquationOfStateOfSolids})(
    y;
    search_interval = (eps(), 2),
    maxiter = 40,
    verbose = false,
)
    v0 = sum(extrema(search_interval)) / 2 * getparam(eos⁻¹.eos).v0  # v0 can be negative
    @assert _ispositive(minimum(v0))  # No negative volume
    v = find_zero(
        guess -> eos⁻¹.eos(guess) - y,
        v0,
        Order16();
        maxevals = maxiter,
        verbose = verbose,
    )
    if !_ispositive(v)
        @warn "the volume found is negative!"
    end
    return v
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
           Base.Fix1(filter, x -> isapprox(eos(x), y; rtol = rtol)) |>
           v -> MVector(v...)
end

# Idea from https://discourse.julialang.org/t/functional-inverse/10959/6
Base.literal_pow(::typeof(^), eos::EquationOfStateOfSolids, ::Val{-1}) = Inverted(eos)

end
