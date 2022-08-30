using Functors: fmap

using .FiniteStrains: EulerianStrainFromVolume, NaturalStrainFromVolume, Eulerian, Natural

import Unitful: ustrip
import .FiniteStrains: straintype

export orderof

include("types.jl")
include("ev.jl")
include("pv.jl")
include("bv.jl")

# Ref: https://github.com/JuliaLang/julia/blob/4a2830a/base/array.jl#L125
"""
    orderof(x::FiniteStrainParameters)

Return the order of a `FiniteStrainParameters`.

# Examples
```jldoctest
julia> orderof(BirchMurnaghan(40, 0.5, 4, 0)) == 3
true
```
"""
orderof(::Type{<:FiniteStrainParameters{N}}) where {N} = N
orderof(x::FiniteStrainParameters) = orderof(typeof(x))

function Base.show(io::IO, params::Parameters)  # See https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(params)
        print(io, typeof(params), "(  ")
        for f in propertynames(params)
            print(io, f, " = ", getproperty(params, f), "  ")
        end
        print(io, ')')
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        println(io, typeof(params))
        for f in propertynames(params)
            println(io, " ", f, " = ", getproperty(params, f))
        end
    end
end
function Base.show(io::IO, eos::EquationOfStateOfSolids)
    params = eos.param
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(eos)
        print(io, typeof(eos), "(  ")
        for f in propertynames(params)
            print(io, f, " = ", getproperty(params, f), "  ")
        end
        print(io, ')')
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        println(io, typeof(eos))
        for f in propertynames(params)
            println(io, " ", f, " = ", getproperty(params, f))
        end
    end
end

straintype(::Type{<:BirchMurnaghan}) = EulerianStrain
straintype(::Type{<:PoirierTarantola}) = NaturalStrain
straintype(x::FiniteStrainParameters) = straintype(typeof(x))

Base.eltype(::Type{<:Parameters{T}}) where {T} = T

"""
    float(p::Parameters)

Convert all elements of a `Parameters` to floating point data types.
"""
Base.float(p::Parameters) = fmap(float, p)  # Not used but may be useful

"""
    isreal(p::Parameters)

Test whether all `p`'s elements are numerically equal to some real number.
"""
Base.isreal(p::Parameters) = all(isreal(getfield(p, i)) for i in 1:nfields(p))  # Not used but may be useful

"""
    real(p::Parameters)

Construct a real `Parameters` from the real parts of the elements of p.
"""
Base.real(p::Parameters) = fmap(real, p)  # Not used but may be useful

"""
    ustrip(p::Parameters)

Strip units from a `Parameters`.
"""
ustrip(p::Parameters) = fmap(ustrip, p)
