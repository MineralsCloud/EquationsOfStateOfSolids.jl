module Collections

using AutoHashEquals: @auto_hash_equals
using UnPack: @unpack

export Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    Eulerian,
    Lagrangian,
    Natural,
    Infinitesimal,
    EnergyEOS,
    PressureEOS,
    BulkModulusEOS,
    orderof,
    nextorder,
    strain_from_volume,
    volume_from_strain

abstract type Parameters{T} end
abstract type FiniteStrainParameters{N,T} <: Parameters{T} end
@auto_hash_equals struct Murnaghan{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct BirchMurnaghan2nd{T} <: FiniteStrainParameters{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
@auto_hash_equals struct BirchMurnaghan3rd{T} <: FiniteStrainParameters{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct BirchMurnaghan4th{T} <: FiniteStrainParameters{4,T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
    BirchMurnaghan4th{T}(v0, b0, b′0, b′′0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b′′0, e0)
end
@auto_hash_equals struct PoirierTarantola2nd{T} <: FiniteStrainParameters{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
@auto_hash_equals struct PoirierTarantola3rd{T} <: FiniteStrainParameters{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct PoirierTarantola4th{T} <: FiniteStrainParameters{4,T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
    PoirierTarantola4th{T}(v0, b0, b′0, b′′0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b′′0, e0)
end
@auto_hash_equals struct Vinet{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Vinet{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
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
    f = strain_from_volume(Eulerian(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * (1 + f * (b′0 - 4))
end
function (eos::EnergyEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return e0 + 3b0 * v0 / 8 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12f * (b′0 - 4) + 12)
end
function (eos::EnergyEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * ((b′0 - 2) * f + 1)
end
function (eos::EnergyEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return e0 + 9b0 * v0 * f^2 * (3f^2 * (h + 3b′0 + 3) + 4f * (b′0 + 2) + 4)
end
function (eos::EnergyEOS{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = 1 - (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end

function (eos::PressureEOS{<:Murnaghan})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function (eos::PressureEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return 3b0 * f * (2f + 1)^(5 / 2)
end
function (eos::PressureEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return 3f / 2 * b0 * (2f + 1)^(5 / 2) * (2 + 3f * (b′0 - 4))
end
function (eos::PressureEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63b′0 + 143) * f^2 + 9f * (b′0 - 4) + 6)
end
function (eos::PressureEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -3b0 * f * exp(-3f)
end
function (eos::PressureEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -3b0 / 2 * f * exp(-3f) * (3f * (b′0 - 2) + 1)
end
function (eos::PressureEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return -3b0 / 2 * f * exp(-3f) * (3f^2 * (h + 3b′0 + 3) + 3f * (b′0 - 2) + 2)
end
function (eos::PressureEOS{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end

(eos::BulkModulusEOS{<:Murnaghan})(v) = eos.param.b0 + PressureEOS(eos.param)(v)
function (eos::BulkModulusEOS{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function (eos::BulkModulusEOS{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (eos::BulkModulusEOS{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return b0 / 6 *
           (2f + 1)^(5 / 2) *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function (eos::BulkModulusEOS{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return b0 * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEOS{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -b0 / 2 * exp(-3f) * (9f^2 * (b′0 - 2) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b′′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    h = b′′0 * b0 + b′0^2
    return -b0 / 2 *
           exp(-3f) *
           (9f^3 * (h + 3b′0 + 3) - 9f^2 * (h + 2b′0 + 1) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEOS{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, ξ = (v / v0)^(1 / 3), 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end

orderof(::FiniteStrainParameters{N}) where {N} = N

abstract type FiniteStrain end  # Trait
struct Eulerian <: FiniteStrain end
struct Lagrangian <: FiniteStrain end
struct Natural <: FiniteStrain end
struct Infinitesimal <: FiniteStrain end

strain_from_volume(::Eulerian, v0) = v -> ((v0 / v)^(2 / 3) - 1) / 2
strain_from_volume(::Lagrangian, v0) = v -> ((v / v0)^(2 / 3) - 1) / 2
strain_from_volume(::Natural, v0) = v -> log(v / v0) / 3
strain_from_volume(::Infinitesimal, v0) = v -> 1 - (v0 / v)^(1 / 3)

volume_from_strain(::Eulerian, v0) = f -> v0 / (2f + 1)^(3 / 2)
volume_from_strain(::Lagrangian, v0) = f -> v0 * (2f + 1)^(3 / 2)
volume_from_strain(::Natural, v0) = f -> v0 * exp(3f)
volume_from_strain(::Infinitesimal, v0) = f -> v0 / (1 - f)^3

function Base.show(io::IO, eos::Parameters)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
    if get(io, :compact, false)
        Base.show_default(IOContext(io, :limit => true), eos)
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        T = typeof(eos)
        println(io, T)
        for f in propertynames(eos)
            println(io, " ", f, " = ", getproperty(eos, f))
        end
    end
end # function Base.show

(::Type{T})(arr::AbstractVector) where {T<:Parameters} = T{eltype(arr)}(arr...)
(::Type{T})(args...) where {T<:Parameters} = T([args...])

end
