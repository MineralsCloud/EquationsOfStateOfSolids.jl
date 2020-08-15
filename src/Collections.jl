module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: AbstractQuantity, @u_str

export BirchMurnaghan,
    Eulerian,
    Lagrangian,
    Natural,
    Infinitesimal,
    energyeos,
    pressureeos,
    bulkmoduluseos,
    nextorder,
    strain_from_volume,
    volume_from_strain

abstract type EossParameters{T} end

abstract type FiniteStrainEossParameters{N,T} <: EossParameters{T} end

struct BirchMurnaghan{N,T} <: FiniteStrainEossParameters{N,T}
    x0::NTuple{N,T}
    e0::T
    function BirchMurnaghan{N,T}(x0, e0) where {N,T}
        @assert 2 <= N <= 5
        new(x0, e0)
    end
end
BirchMurnaghan(x0::NTuple{N,T}, e0 = zero(x0[1] * x0[2])) where {N,T} =
    BirchMurnaghan{N,T}(x0, e0)
BirchMurnaghan(arr::AbstractArray{T}) where {T} =
    BirchMurnaghan(Tuple(arr[1:end-1]), arr[end])
BirchMurnaghan(args...) = BirchMurnaghan([args...])

abstract type EquationOfStateOfSolids{T<:EossParameters} end
struct EnergyEoss{T} <: EquationOfStateOfSolids{T}
    params::T
end
struct PressureEoss{T} <: EquationOfStateOfSolids{T}
    params::T
end
struct BulkModulusEoss{T} <: EquationOfStateOfSolids{T}
    params::T
end

energyeos(p::EossParameters) = EnergyEoss(p)
pressureeos(p::EossParameters) = PressureEoss(p)
bulkmoduluseos(p::EossParameters) = BulkModulusEoss(p)

function (eos::EnergyEoss{<:BirchMurnaghan{3}})(v, e0 = zero(eos.params.x0[1] * eos.params.x0[3]))
    v0, b0, b′0 = eos.params.x0
    x = cbrt(v0 / v)
    y = x^2 - 1
    return 9 / 16 * b0 * v0 * y^2 * (6 - 4 * x^2 + b′0 * y) + e0
end

function (eos::PressureEoss{<:BirchMurnaghan{3}})(v)
    v0, b0, b′0 = eos.params.x0
    f = strain_from_volume(Eulerian(), v0)(v)
    return 3f / 2 * b0 * sqrt(2f + 1)^5 * (2 + 3f * (b′0 - 4))
end

nextorder(::Type{BirchMurnaghan{N}}) where {N} = BirchMurnaghan{N + 1}

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

fieldvalues(x::EossParameters) = x.x0

function Base.getproperty(value::FiniteStrainEossParameters, name::Symbol)
    if name == :v0
        return value.x0[1]
    elseif name == :b0
        return value.x0[2]
    elseif name == :b′0
        return value.x0[3]
    else
        return getfield(value, name)
    end
end

end
