module Collections

using AutoHashEquals: @auto_hash_equals
using Unitful: AbstractQuantity, NoUnits, ħ, me
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
    EnergyEquation,
    PressureEquation,
    BulkModulusEquation,
    orderof,
    atomic_number,
    volume2strain,
    strain2volume,
    straintype,
    getparam

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
julia> orderof(BirchMurnaghan3rd(40, 0.5, 4, 0)) == 3
true
```
"""
orderof(::Type{<:FiniteStrainParameters{N}}) where {N} = N
orderof(x::FiniteStrainParameters) = orderof(typeof(x))

atomic_number(::Type{<:Holzapfel{Z}}) where {Z} = Z
atomic_number(x::Holzapfel) = atomic_number(typeof(x))

abstract type FiniteStrain end  # Trait
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

"""
    volume2strain(::EulerianStrain, v0)
    volume2strain(::LagrangianStrain, v0)
    volume2strain(::NaturalStrain, v0)
    volume2strain(::InfinitesimalStrain, v0)

Return a function of `v` that calculates the `FiniteStrain` from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
volume2strain(::EulerianStrain, v0) = v -> ((v0 / v)^_⅔ - 1) / 2
volume2strain(::LagrangianStrain, v0) = v -> ((v / v0)^_⅔ - 1) / 2
volume2strain(::NaturalStrain, v0) = v -> log(v / v0) / 3
volume2strain(::InfinitesimalStrain, v0) = v -> 1 - (v0 / v)^_⅓

"""
    strain2volume(::EulerianStrain, v0)
    strain2volume(::LagrangianStrain, v0)
    strain2volume(::NaturalStrain, v0)
    strain2volume(::InfinitesimalStrain, v0)

Return a function of `f` that calculates the corresponding volume from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
function strain2volume(::EulerianStrain, v0)
    return function (f)
        if isreal(f)
            f = real(f)
        else
            throw(DomainError("strain $f is complex!"))
        end
        if f < -1 / 2
            throw(DomainError("strain $f < -0.5! Volume will be complex!"))
        else
            return v0 / (2f + 1)^_1½
        end
    end
end
function strain2volume(::LagrangianStrain, v0)
    return function (f)
        if isreal(f)
            f = real(f)
        else
            throw(DomainError("strain $f is complex!"))
        end
        if f < -1 / 2
            throw(DomainError("strain $f < -0.5! Volume will be complex!"))
        else
            return v0 * (2f + 1)^_1½
        end
    end
end
strain2volume(::NaturalStrain, v0) = f -> v0 * exp(3f)
strain2volume(::InfinitesimalStrain, v0) = f -> v0 / (1 - f)^3

"""
    Dⁿᵥf(s::EulerianStrain, deg, v0)
    Dⁿᵥf(s::LagrangianStrain, deg, v0)
    Dⁿᵥf(s::NaturalStrain, deg, v0)
    Dⁿᵥf(s::InfinitesimalStrain, deg, v0)

Return a function of `v` that calculates the `deg`th order derivative of strain wrt volume from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
function Dⁿᵥf(s::EulerianStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v0 / v)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 1) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::LagrangianStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v / v0)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 5) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::NaturalStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return 1 / 3 / v
        else  # Recursion
            return -(deg - 1) / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::InfinitesimalStrain, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return (1 - volume2strain(s, v0)(v))^4 / 3 / v0
        else  # Recursion
            return -(3deg - 2) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end

straintype(::Type{<:BirchMurnaghan}) = EulerianStrain
straintype(::Type{<:PoirierTarantola}) = NaturalStrain
straintype(x::FiniteStrainParameters) = straintype(typeof(x))

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

(::Type{T})(arr::AbstractVector) where {T<:Parameters} = T{eltype(arr)}(arr...)
(::Type{T})(args...) where {T<:Parameters} = T([args...])
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

Base.eltype(::Type{<:Parameters{T}}) where {T} = T

end
