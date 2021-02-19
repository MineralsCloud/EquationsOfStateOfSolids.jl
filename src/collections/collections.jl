using AutoHashEquals: @auto_hash_equals
using ConstructionBase: constructorof
using UnPack: @unpack

using .FiniteStrains: ToEulerianStrain, ToNaturalStrain, EulerianStrain, NaturalStrain

import .FiniteStrains: straintype

export Murnaghan1st,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    EnergyEquation,
    PressureEquation,
    BulkModulusEquation,
    EnergyFrom,
    PressureFrom,
    BulkModulusFrom,
    orderof,
    getparam

include("types.jl")
include("ev.jl")
include("pv.jl")
include("bv.jl")

_ispositive(x) = x > zero(x)  # Do not export!

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

function (::Type{T})(args...) where {T<:Parameters}
    E = Base.promote_typeof(args...)
    return constructorof(T){E}(args...)  # Cannot use `T.(args...)`! For `AbstractQuantity` they will fail!
end
function Murnaghan(args...)
    N = length(args)
    if N == 4
        return Murnaghan1st(args...)
    elseif N == 5
        return Murnaghan2nd(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end
"""
    BirchMurnaghan(args...)

Construct a `BirchMurnaghan` based on the length of arguments, where `e0` must be provided.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan3rd`](@ref), [`BirchMurnaghan4th`](@ref)
"""
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
"""
    PoirierTarantola(args...)

Construct a `PoirierTarantola` based on the length of arguments, where `e0` must be provided.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref), [`PoirierTarantola4th`](@ref)
"""
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

straintype(::Type{<:BirchMurnaghan}) = EulerianStrain
straintype(::Type{<:PoirierTarantola}) = NaturalStrain
straintype(x::FiniteStrainParameters) = straintype(typeof(x))

Base.eltype(::Type{<:Parameters{T}}) where {T} = T
