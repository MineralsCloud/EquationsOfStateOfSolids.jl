module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: AbstractQuantity, NoUnits, ħ, me
using UnPack: @unpack

export Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    AntonSchmidt,
    Holzapfel,
    EulerianStrain,
    LagrangianStrain,
    NaturalStrain,
    InfinitesimalStrain,
    EnergyEOS,
    PressureEOS,
    BulkModulusEOS,
    orderof,
    atomic_number,
    volume2strain,
    strain2volume,
    straintype,
    getparam

const FERMI_GAS_CONSTANT = (3π^2)^(2 / 3) * ħ^2 / 5 / me

abstract type Power end
struct TwoThirds <: Power end
struct OneThird <: Power end
struct FiveHalves <: Power end
struct ThreeHalves <: Power end

const _⅔ = TwoThirds()
const _⅓ = OneThird()
const _2½ = FiveHalves()
const _1½ = ThreeHalves()

Base.:(^)(x, ::TwoThirds) = x^(2 // 3)
Base.:(^)(x::Union{Real,AbstractQuantity}, ::TwoThirds) = cbrt(x^2)
Base.:(^)(x, ::OneThird) = x^(1 // 3)
Base.:(^)(x::Union{Real,AbstractQuantity}, ::OneThird) = cbrt(x)
Base.:(^)(x, ::FiveHalves) = x^(5 // 2)
Base.:(^)(x::Union{Real,AbstractQuantity}, ::FiveHalves) = sqrt(x^5)
Base.:(^)(x, ::ThreeHalves) = x^(3 // 2)
Base.:(^)(x::Union{Real,AbstractQuantity}, ::ThreeHalves) = sqrt(x^3)

abstract type Parameters{T} end
abstract type FiniteStrainParameters{N,T} <: Parameters{T} end
abstract type BirchMurnaghan{N,T} <: FiniteStrainParameters{N,T} end
abstract type PoirierTarantola{N,T} <: FiniteStrainParameters{N,T} end
"""
    Murnaghan(v0, b0, b′0, e0=zero(v0 * b0))

Create a Murnaghan equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct Murnaghan{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    BirchMurnaghan2nd(v0, b0, e0=zero(v0 * b0))

Create a Birch–Murnaghan second order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct BirchMurnaghan2nd{T} <: BirchMurnaghan{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
"""
    BirchMurnaghan3rd(v0, b0, b′0, e0=zero(v0 * b0))

Create a Birch–Murnaghan third order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct BirchMurnaghan3rd{T} <: BirchMurnaghan{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    BirchMurnaghan4th(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Birch–Murnaghan fourth order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure. "`″`" can be typed by `\\pprime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct BirchMurnaghan4th{T} <: BirchMurnaghan{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    BirchMurnaghan4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
"""
    PoirierTarantola2nd(v0, b0, e0=zero(v0 * b0))

Create a Poirier–Tarantola second order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola2nd{T} <: PoirierTarantola{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
"""
    PoirierTarantola3rd(v0, b0, b′0, e0=zero(v0 * b0))

Create a Poirier–Tarantola third order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola3rd{T} <: PoirierTarantola{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    PoirierTarantola4th(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Poirier–Tarantola fourth order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure. "`″`" can be typed by `\\pprime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola4th{T} <: PoirierTarantola{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    PoirierTarantola4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
"""
    Vinet(v0, b0, b′0, e0=zero(v0 * b0))

Create a Vinet equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct Vinet{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Vinet{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct AntonSchmidt{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e∞::T
    AntonSchmidt{T}(v0, b0, b′0, e∞ = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e∞)
end
@auto_hash_equals struct Holzapfel{Z,T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    function Holzapfel{Z,T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {Z,T}
        @assert 1 <= Z <= 118 "elements are between 1 and 118!"
        return new(v0, b0, b′0, e0)
    end
end

abstract type EquationOfState{T<:Parameters} end
abstract type EquationOfStateOfSolids{T} <: EquationOfState{T} end
struct EnergyEOS{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct PressureEOS{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct BulkModulusEOS{T} <: EquationOfStateOfSolids{T}
    param::T
end

function (eos::EnergyEOS{<:Murnaghan})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    x, y = b′0 - 1, (v0 / v)^b′0
    return e0 + b0 / b′0 * v * (y / x + 1) - v0 * b0 / x
end
function (eos::EnergyEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0, e0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * (1 + f * (b′0 - 4))
end
function (eos::EnergyEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 3b0 * v0 / 8 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12f * (b′0 - 4) + 12)
end
function (eos::EnergyEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
"""
    (eos::EnergyEOS{<:PoirierTarantola3rd})(v)

Evaluate a Poirier--Tarantola energy EOS on volume `v`.

```math
E(V) = E_0 + \\frac {B_0 V_0}{2}\\left[ \\ln \\left( \\frac{V_0}{V}\\right) \\right]^{2} + \\frac{B_0 V_0}{6} \\left[ \\ln \\left(\\frac{V_0}V{} \\right) \\right]^3 \\left(B_0^\\prime - 2 \\right)
```
"""
function (eos::EnergyEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * ((2 - b′0) * f + 1)
end
function (eos::EnergyEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 9b0 * v0 * f^2 * (3f^2 * (h + 3b′0 + 3) + 4f * (b′0 + 2) + 4)
end
function (eos::EnergyEOS{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    x, y = 1 - (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end
function (f::EnergyEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0, e∞ = getparam(eos)
    n₊₁ = -b′0 / 2 + 1
    x = v / v0
    return e∞ + b0 * v0 / n₊₁ * x^n₊₁ * (log(x) - 1 / n₊₁)
end

function (eos::PressureEOS{<:Murnaghan})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function (eos::PressureEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return 3b0 * f * (2f + 1)^_2½
end
function (eos::PressureEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return 3f / 2 * b0 * (2f + 1)^_2½ * (2 + 3f * (b′0 - 4))
end
function (eos::PressureEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 2 * (2f + 1)^_2½ * ((9h - 63b′0 + 143) * f^2 + 9f * (b′0 - 4) + 6)
end
function (eos::PressureEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return -3b0 * f * exp(-3f)
end
"""
    (eos::PressureEOS{<:PoirierTarantola3rd})(v)

Evaluate a Poirier--Tarantola pressure EOS on volume `v`.

```math
P(V) = B_0 \\frac{V_0}{V} \\left[\\ln \\left( \\frac{V_0}{V} \\right) + \\frac{\\left( B_0^\\prime -2 \\right) }{2} \\left[ \\ln \\left( \\frac{V_0}{V} \\right) \\right]^2\\right]
```
"""
function (eos::PressureEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return -3b0 / 2 * f * exp(-3f) * (3f * (b′0 - 2) + 1)
end
function (eos::PressureEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return -3b0 / 2 * f * exp(-3f) * (3f^2 * (h + 3b′0 + 3) + 3f * (b′0 - 2) + 2)
end
function (eos::PressureEOS{<:Vinet})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, y = (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end
function (f::PressureEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, n = v / v0, -b′0 / 2
    return -b0 * x^n * log(x)
end
function (eos::PressureEOS{<:Holzapfel{Z}})(v) where {Z}
    @unpack v0, b0, b′0 = getparam(eos)
    η = (v / v0)^_⅓
    p0 = FERMI_GAS_CONSTANT * (Z / v0)^(5 / 3)
    c0 = -log(3b0 / p0 |> NoUnits)
    c2 = 3 / 2 * (b′0 - 3) - c0
    return 3b0 * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
end

(eos::BulkModulusEOS{<:Murnaghan})(v) = getparam(eos).b0 + PressureEOS(getparam(eos))(v)
function (eos::BulkModulusEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return b0 * (7f + 1) * (2f + 1)^_2½
end
function (eos::BulkModulusEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    return b0 / 2 * (2f + 1)^_2½ * ((27f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (eos::BulkModulusEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 6 *
           (2f + 1)^_2½ *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function (eos::BulkModulusEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return b0 * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    return -b0 / 2 * exp(-3f) * (9f^2 * (b′0 - 2) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return -b0 / 2 *
           exp(-3f) *
           (9f^3 * (h + 3b′0 + 3) - 9f^2 * (h + 2b′0 + 1) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:Vinet})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, ξ = (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end
function (f::BulkModulusEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, n = v / v0, -b′0 / 2
    return b0 * x^n * (1 + n * log(x))
end

# Ref: https://github.com/JuliaLang/julia/blob/4a2830a/base/array.jl#L125
"""
    orderof(x::FiniteStrainParameters)

Return the order of a `FiniteStrainParameters`.

# Examples
```jldoctest
julia> orderof(BirchMurnaghan3rd(40, 0.5, 4, 0)) == 3
true
```
"""
orderof(::Type{<:FiniteStrainParameters{N}}) where {N} = N
orderof(x::FiniteStrainParameters) = orderof(typeof(x))

atomic_number(::Type{<:Holzapfel{Z}}) where {Z} = Z
atomic_number(x::Holzapfel) = atomic_number(typeof(x))

abstract type FiniteStrain end  # Trait
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

"""
    volume2strain(::EulerianStrain, v0)
    volume2strain(::LagrangianStrain, v0)
    volume2strain(::NaturalStrain, v0)
    volume2strain(::InfinitesimalStrain, v0)

Return a function of `v` that calculates the `FiniteStrain` from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
volume2strain(::EulerianStrain, v0) = v -> ((v0 / v)^_⅔ - 1) / 2
volume2strain(::LagrangianStrain, v0) = v -> ((v / v0)^_⅔ - 1) / 2
volume2strain(::NaturalStrain, v0) = v -> log(v / v0) / 3
volume2strain(::InfinitesimalStrain, v0) = v -> 1 - (v0 / v)^_⅓

"""
    strain2volume(::EulerianStrain, v0)
    strain2volume(::LagrangianStrain, v0)
    strain2volume(::NaturalStrain, v0)
    strain2volume(::InfinitesimalStrain, v0)

Return a function of `f` that calculates the corresponding volume from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
strain2volume(::EulerianStrain, v0) = f -> v0 / (2f + 1)^_1½
strain2volume(::LagrangianStrain, v0) = f -> v0 * (2f + 1)^_1½
strain2volume(::NaturalStrain, v0) = f -> v0 * exp(3f)
strain2volume(::InfinitesimalStrain, v0) = f -> v0 / (1 - f)^3

"""
    Dⁿᵥf(s::EulerianStrain, deg, v0)
    Dⁿᵥf(s::LagrangianStrain, deg, v0)
    Dⁿᵥf(s::NaturalStrain, deg, v0)
    Dⁿᵥf(s::InfinitesimalStrain, deg, v0)

Return a function of `v` that calculates the `deg`th order derivative of strain wrt volume from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
function Dⁿᵥf(s::EulerianStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v0 / v)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 1) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::LagrangianStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v / v0)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 5) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::NaturalStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return 1 / 3 / v
        else  # Recursion
            return -(deg - 1) / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::InfinitesimalStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return (1 - volume2strain(s, v0)(v))^4 / 3 / v0
        else  # Recursion
            return -(3deg - 2) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end

straintype(::Type{<:BirchMurnaghan}) = EulerianStrain
straintype(::Type{<:PoirierTarantola}) = NaturalStrain
straintype(x::FiniteStrainParameters) = straintype(typeof(x))

function Base.show(io::IO, param::Parameters)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
    if get(io, :compact, false)
        Base.show_default(IOContext(io, :limit => true), param)
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        T = typeof(param)
        println(io, T)
        for f in propertynames(param)
            println(io, " ", f, " = ", getproperty(param, f))
        end
    end
end # function Base.show

(::Type{T})(arr::AbstractVector) where {T<:Parameters} = T{eltype(arr)}(arr...)
(::Type{T})(args...) where {T<:Parameters} = T([args...])
function BirchMurnaghan(args...)
    N = length(args)
    if N == 3
        return BirchMurnaghan2nd(args...)
    elseif N == 4
        return BirchMurnaghan3rd(args...)
    elseif N == 5
        return BirchMurnaghan4th(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end
function PoirierTarantola(args...)
    N = length(args)
    if N == 3
        return PoirierTarantola2nd(args...)
    elseif N == 4
        return PoirierTarantola3rd(args...)
    elseif N == 5
        return PoirierTarantola4th(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end

"""
    getparam(eos::EquationOfStateOfSolids)

Get the `Parameters` from an `EquationOfStateOfSolids`.
"""
getparam(eos::EquationOfStateOfSolids) = eos.param

Base.eltype(::Type{<:Parameters{T}}) where {T} = T

end
