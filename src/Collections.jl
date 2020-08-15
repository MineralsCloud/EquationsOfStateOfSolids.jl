module Collections

using AutoHashEquals: @auto_hash_equals
using UnPack: @unpack

export BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    Vinet,
    Eulerian,
    Lagrangian,
    Natural,
    Infinitesimal,
    EnergyEoss,
    PressureEoss,
    BulkModulusEoss,
    orderof,
    nextorder,
    strain_from_volume,
    volume_from_strain

abstract type EossParam{T} end
abstract type FiniteStrainEossParam{N,T} <: EossParam{T} end
struct BirchMurnaghan2nd{T} <: FiniteStrainEossParam{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
struct BirchMurnaghan3rd{T} <: FiniteStrainEossParam{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
struct PoirierTarantola2nd{T} <: FiniteStrainEossParam{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
struct PoirierTarantola3rd{T} <: FiniteStrainEossParam{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
struct Vinet{T} <: EossParam{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Vinet{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end

abstract type EquationOfStateOfSolids{T<:EossParam} end
struct EnergyEoss{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct PressureEoss{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct BulkModulusEoss{T} <: EquationOfStateOfSolids{T}
    param::T
end

function (eos::EnergyEoss{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEoss{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * (1 + f * (b′0 - 4))
end
function (eos::EnergyEoss{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEoss{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * ((b′0 - 2) * f + 1)
end
function (eos::EnergyEoss{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = 1 - cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end

function (eos::PressureEoss{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return 3b0 * f * (2f + 1)^(5 / 2)
end
function (eos::PressureEoss{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return 3f / 2 * b0 * sqrt(2f + 1)^5 * (2 + 3f * (b′0 - 4))
end
function (eos::PressureEoss{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -3b0 * f * exp(-3f)
end
function (eos::PressureEoss{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -3b0 / 2 * f * exp(-3f) * (3f * (b′0 - 2) + 1)
end
function (eos::PressureEoss{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, y = cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end

function (eos::BulkModulusEoss{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function (eos::BulkModulusEoss{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Eulerian(), v0)(v)
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (eos::BulkModulusEoss{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return b0 * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEoss{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    f = strain_from_volume(Natural(), v0)(v)
    return -b0 / 2 * exp(-3f) * (9f^2 * (b′0 - 2) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEoss{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = eos.param
    x, ξ = cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end

orderof(::FiniteStrainEossParam{N}) where {N} = N

abstract type FiniteStrain end  # Trait
struct Eulerian <: FiniteStrain end
struct Lagrangian <: FiniteStrain end
struct Natural <: FiniteStrain end
struct Infinitesimal <: FiniteStrain end

strain_from_volume(::Eulerian, v0) = v -> (cbrt(v0 / v)^2 - 1) / 2
strain_from_volume(::Lagrangian, v0) = v -> (cbrt(v / v0)^2 - 1) / 2
strain_from_volume(::Natural, v0) = v -> log(v / v0) / 3
strain_from_volume(::Infinitesimal, v0) = v -> 1 - cbrt(v0 / v)

volume_from_strain(::Eulerian, v0) = f -> v0 / (2f + 1)^(3 / 2)
volume_from_strain(::Lagrangian, v0) = f -> v0 * (2f + 1)^(3 / 2)
volume_from_strain(::Natural, v0) = f -> v0 * exp(3f)
volume_from_strain(::Infinitesimal, v0) = f -> v0 / (1 - f)^3

function Base.show(io::IO, eos::EossParam)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
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

for T in (
    :BirchMurnaghan2nd,
    :BirchMurnaghan3rd,
    :PoirierTarantola2nd,
    :PoirierTarantola3rd,
    :Vinet,
)
    eval(quote
        $T(arr::AbstractVector) = $T{eltype(arr)}(arr...)
        $T(args...) = $T([args...])
    end)
end

end
