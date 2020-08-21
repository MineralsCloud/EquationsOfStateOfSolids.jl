module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: NoUnits, ħ, me
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
    whatstrain

const FERMI_GAS_CONSTANT = (3π^2)^(2 / 3) * ħ^2 / 5 / me

abstract type Parameters{T} end
abstract type FiniteStrainParameters{N,T} <: Parameters{T} end
abstract type BirchMurnaghan{N,T} <: FiniteStrainParameters{N,T} end
abstract type PoirierTarantola{N,T} <: FiniteStrainParameters{N,T} end
@auto_hash_equals struct Murnaghan{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct BirchMurnaghan2nd{T} <: BirchMurnaghan{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
@auto_hash_equals struct BirchMurnaghan3rd{T} <: BirchMurnaghan{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct BirchMurnaghan4th{T} <: BirchMurnaghan{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    BirchMurnaghan4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
@auto_hash_equals struct PoirierTarantola2nd{T} <: PoirierTarantola{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
@auto_hash_equals struct PoirierTarantola3rd{T} <: PoirierTarantola{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct PoirierTarantola4th{T} <: PoirierTarantola{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    PoirierTarantola4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
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
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = b′0 - 1, (v0 / v)^b′0
    return e0 + b0 / b′0 * v * (y / x + 1) - v0 * b0 / x
end
function (eos::EnergyEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * (1 + f * (b′0 - 4))
end
function (eos::EnergyEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 3b0 * v0 / 8 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12f * (b′0 - 4) + 12)
end
function (eos::EnergyEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * ((b′0 - 2) * f + 1)
end
function (eos::EnergyEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 9b0 * v0 * f^2 * (3f^2 * (h + 3b′0 + 3) + 4f * (b′0 + 2) + 4)
end
function (eos::EnergyEOS{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = 1 - (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end
function (f::EnergyEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0, e∞ = eos.param
    n₊₁ = -b′0 / 2 + 1
    x = v / v0
    return e∞ + b0 * v0 / n₊₁ * x^n₊₁ * (log(x) - 1 / n₊₁)
end

function (eos::PressureEOS{<:Murnaghan})(v)
    @unpack v0, b0, b′0 = eos.param
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function (eos::PressureEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return 3b0 * f * (2f + 1)^(5 / 2)
end
function (eos::PressureEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return 3f / 2 * b0 * (2f + 1)^(5 / 2) * (2 + 3f * (b′0 - 4))
end
function (eos::PressureEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63b′0 + 143) * f^2 + 9f * (b′0 - 4) + 6)
end
function (eos::PressureEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return -3b0 * f * exp(-3f)
end
function (eos::PressureEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return -3b0 / 2 * f * exp(-3f) * (3f * (b′0 - 2) + 1)
end
function (eos::PressureEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return -3b0 / 2 * f * exp(-3f) * (3f^2 * (h + 3b′0 + 3) + 3f * (b′0 - 2) + 2)
end
function (eos::PressureEOS{<:Vinet})(v)
    @unpack v0, b0, b′0 = eos.param
    x, y = (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end
function (f::PressureEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = eos.param
    x, n = v / v0, -b′0 / 2
    return -b0 * x^n * log(x)
end
function (eos::PressureEOS{<:Holzapfel{Z}})(v) where {Z}
    @unpack v0, b0, b′0 = eos.param
    η = (v / v0)^(1 / 3)
    p0 = FERMI_GAS_CONSTANT * (Z / v0)^(5 / 3)
    c0 = -log(3b0 / p0 |> NoUnits)
    c2 = 3 / 2 * (b′0 - 3) - c0
    return 3b0 * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
end

(eos::BulkModulusEOS{<:Murnaghan})(v) = eos.param.b0 + PressureEOS(eos.param)(v)
function (eos::BulkModulusEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function (eos::BulkModulusEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (eos::BulkModulusEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = eos.param
    f = volume2strain(EulerianStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 6 *
           (2f + 1)^(5 / 2) *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function (eos::BulkModulusEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return b0 * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    return -b0 / 2 * exp(-3f) * (9f^2 * (b′0 - 2) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = eos.param
    f = volume2strain(NaturalStrain(), v0)(v)
    h = b″0 * b0 + b′0^2
    return -b0 / 2 *
           exp(-3f) *
           (9f^3 * (h + 3b′0 + 3) - 9f^2 * (h + 2b′0 + 1) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:Vinet})(v)
    @unpack v0, b0, b′0 = eos.param
    x, ξ = (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end
function (f::BulkModulusEOS{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = eos.param
    x, n = v / v0, -b′0 / 2
    return b0 * x^n * (1 + n * log(x))
end

# Ref: https://github.com/JuliaLang/julia/blob/4a2830a/base/array.jl#L125
orderof(::Type{<:FiniteStrainParameters{N}}) where {N} = N
orderof(x::FiniteStrainParameters) = orderof(typeof(x))

atomic_number(::Type{<:Holzapfel{Z}}) where {Z} = Z
atomic_number(x::Holzapfel) = atomic_number(typeof(x))

abstract type FiniteStrain end  # Trait
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

volume2strain(::EulerianStrain, v0) = v -> ((v0 / v)^(2 / 3) - 1) / 2
volume2strain(::LagrangianStrain, v0) = v -> ((v / v0)^(2 / 3) - 1) / 2
volume2strain(::NaturalStrain, v0) = v -> log(v / v0) / 3
volume2strain(::InfinitesimalStrain, v0) = v -> 1 - (v0 / v)^(1 / 3)

strain2volume(::EulerianStrain, v0) = f -> v0 / (2f + 1)^(3 / 2)
strain2volume(::LagrangianStrain, v0) = f -> v0 * (2f + 1)^(3 / 2)
strain2volume(::NaturalStrain, v0) = f -> v0 * exp(3f)
strain2volume(::InfinitesimalStrain, v0) = f -> v0 / (1 - f)^3

function strain_volume_derivative(s::EulerianStrain, v0, v, deg::Integer)
    if deg == 1
        return -(v0 / v)^(2 / 3) / 3 / v
    else  # Recursion
        return -(3deg - 1) / 3 / v * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::LagrangianStrain, v0, v, deg::Integer)
    if deg == 1
        return -(v / v0)^(2 / 3) / 3 / v
    else  # Recursion
        return -(3deg - 5) / 3 / v * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::NaturalStrain, v0, v, deg::Integer)
    if deg == 1
        return 1 / 3 / v
    else  # Recursion
        return -(deg - 1) / v * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::InfinitesimalStrain, v0, v, deg::Integer)
    if deg == 1
        return (1 - volume2strain(s, v0)(v))^4 / 3 / v0
    else  # Recursion
        return -(3deg - 2) / 3 / v * strain_volume_derivative(s, v0, v, deg - 1)
    end
end

whatstrain(::Type{<:BirchMurnaghan}) = EulerianStrain()
whatstrain(::Type{<:PoirierTarantola}) = NaturalStrain()
whatstrain(x::FiniteStrainParameters) = whatstrain(typeof(x))

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

Base.eltype(::Type{<:Parameters{T}}) where {T} = T

end
