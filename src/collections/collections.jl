using Functors: fmap
using UnPack: @unpack

using .FiniteStrains: ToEulerianStrain, ToNaturalStrain, EulerianStrain, NaturalStrain

import Unitful: ustrip
import .FiniteStrains: straintype

export getparam, orderof

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
end

"""
    getparam(eos::EquationOfStateOfSolids)

Get the `Parameters` from an `EquationOfStateOfSolids`.
"""
getparam(eos::EquationOfStateOfSolids) = eos.param

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
