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

const BirchMurnaghan3rd = BirchMurnaghan{3}

function energyeq(p::BirchMurnaghan{3})
    v0, b0, b′0 = p.x0
    function (v)
        f = cbrt(v0 / v)^2 - 1
        return 9f^2 / 2 * v0 * b0 * (1 + (b′0 - 4) * f)
    end
end

function pressureeq(p::BirchMurnaghan{3})
    v0, b0, b′0 = p.x0
    function (v)
        f = cbrt(v0 / v)^2 - 1
        return 3f / 2 * b0 * sqrt(2f + 1)^5 * (2 + 3f * (b′0 - 4))
    end
end

nextorder(::Type{BirchMurnaghan{N}}) where {N} = BirchMurnaghan{N + 1}


end
