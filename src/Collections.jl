module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: AbstractQuantity, @u_str

export BirchMurnaghan, BirchMurnaghan3rd, energyeq, pressureeq, bulkmoduluseq, nextorder

abstract type Parameters{T} end

abstract type ParametersFiniteStrain{N,T} <: Parameters{T} end

struct BirchMurnaghan{N,T} <: ParametersFiniteStrain{N,T}
    x0::NTuple{N,T}
end
function BirchMurnaghan{3}(v0, b0, b′0)
    T = Base.promote_typeof(v0, b0, b′0)
    return BirchMurnaghan{3,T}(Tuple(convert(T, x) for x in (v0, b0, b′0)))
end
# BirchMurnaghan{3}(v0::Real, b0::Real, b′0::Real) = BirchMurnaghan{3}(v0, b0, b′0, 0)
# BirchMurnaghan{3}(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
#     BirchMurnaghan{3}(v0, b0, b′0, 0 * u"eV")

abstract type EquationOfStateOfSolids{T<:Parameters} end
struct EnergyEOSS{T} <: EquationOfStateOfSolids{T}
    params::T
end
struct PressureEOSS{T} <: EquationOfStateOfSolids{T}
    params::T
end
struct BulkModulusEOSS{T} <: EquationOfStateOfSolids{T}
    params::T
end

const BirchMurnaghan3rd = BirchMurnaghan{3}

energyeq(p::Parameters) = EnergyEOSS(p)
pressureeq(p::Parameters) = PressureEOSS(p)
bulkmoduluseq(p::Parameters) = BulkModulusEOSS(p)

function (eos::EnergyEOSS{<:BirchMurnaghan3rd})(v)
    v0, b0, b′0 = eos.params.x0
    x = cbrt(v0 / v)
    y = x^2 - 1
    return 9 / 16 * b0 * v0 * y^2 * (6 - 4 * x^2 + b′0 * y)
end

# function (eos::PressureEOSS{<:BirchMurnaghan3rd})(v)
#     v0, b0, b′0 = eos.params.x0
#     x = cbrt(v0 / v)
#     return 3 / 2 * b0 * (x^7 - x^5) * (1 + 3 / 4 * (b′0 - 4) * (x^2 - 1))
# end

nextorder(::Type{BirchMurnaghan{N}}) where {N} = BirchMurnaghan{N + 1}


end
