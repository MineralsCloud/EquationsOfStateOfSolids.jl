module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: AbstractQuantity, @u_str

export BirchMurnaghan3rd,
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

struct BirchMurnaghan3rd{T} <: FiniteStrainEossParameters{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
BirchMurnaghan3rd(arr::AbstractArray) = BirchMurnaghan3rd{eltype(arr)}(arr...)
BirchMurnaghan3rd(args...) = BirchMurnaghan3rd([args...])

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

function Base.show(io::IO, eos::EossParameters)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
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

end
