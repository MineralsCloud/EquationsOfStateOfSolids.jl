module FiniteStrains

using Unitful: AbstractQuantity, NoUnits

using ..EquationsOfStateOfSolids:
    FiniteStrainParameters, BirchMurnaghan, PoirierTarantola, _⅔, _⅓, _1½

export EulerianStrain,
    LagrangianStrain,
    NaturalStrain,
    InfinitesimalStrain,
    volume2strain,
    strain2volume,
    straintype

abstract type FiniteStrain end  # Trait
# struct EulerianStrain <: FiniteStrain end
# struct LagrangianStrain <: FiniteStrain end
# struct NaturalStrain <: FiniteStrain end
# struct InfinitesimalStrain <: FiniteStrain end
# const Eulerian = EulerianStrain
# const Lagrangian = LagrangianStrain
# const Natural = NaturalStrain
# const Infinitesimal = InfinitesimalStrain

abstract type VolumeToStrain{T} end

struct VolumeToEulerianStrain{T} <: VolumeToStrain{T}
    v0::T
end
(x::VolumeToEulerianStrain)(v) = ((x.v0 / v)^_⅔ - 1) / 2

struct VolumeToLagrangianStrain{T} <: VolumeToStrain{T}
    v0::T
end
(x::VolumeToLagrangianStrain)(v) = ((v / x.v0)^_⅔ - 1) / 2

struct VolumeToNaturalStrain{T} <: VolumeToStrain{T}
    v0::T
end
(x::VolumeToNaturalStrain)(v) = log(v / x.v0) / 3

struct VolumeToInfinitesimalStrain{T} <: VolumeToStrain{T}
    v0::T
end
(::VolumeToInfinitesimalStrain)(v) = 1 - (x.v0 / v)^_⅓

"""
    volume2strain(::EulerianStrain, v0)
    volume2strain(::LagrangianStrain, v0)
    volume2strain(::NaturalStrain, v0)
    volume2strain(::InfinitesimalStrain, v0)

Return a function of `v` that calculates the `FiniteStrain` from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""

"""
    strain2volume(::EulerianStrain, v0)
    strain2volume(::LagrangianStrain, v0)
    strain2volume(::NaturalStrain, v0)
    strain2volume(::InfinitesimalStrain, v0)

Return a function of `f` that calculates the corresponding volume from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
abstract type StrainToVolume{T} end
(x::StrainToVolume)(f::AbstractQuantity) = x(NoUnits(f))

struct EulerianStrainToVolume{T} <: StrainToVolume{T}
    v0::T
end
(x::EulerianStrainToVolume)(f::Real) =
    f < -1 / 2 ? throw(DomainError("strain $f < -0.5! Volume will be complex!")) :
    _EulerianStrainToVolume(x.v0, f)
(x::EulerianStrainToVolume)(f) = _EulerianStrainToVolume(x.v0, f)  # For symbols, etc.
function (x::NaturalStrainToVolume)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(2f + 1)) / 120)  # `(2f + 1)^(3 / 2)` is real for some complex `f`
        return real(_EulerianStrainToVolume(x.v0, f))
    else
        throw(DomainError("volume will be complex!"))
    end
end
_EulerianStrainToVolume(v0, f) = v0 / (2f + 1)^_1½

struct LagrangianStrainToVolume{T} <: StrainToVolume{T}
    v0::T
end
(x::LagrangianStrainToVolume)(f::Real) =
    f < -1 / 2 ? throw(DomainError("strain $f < -0.5! Volume will be complex!")) :
    _LagrangianStrainToVolume(x.v0, f)
(x::LagrangianStrainToVolume)(f) = _LagrangianStrainToVolume(x.v0, f)
function (x::LagrangianStrainToVolume)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(2f + 1)) / 120)  # `(2f + 1)^(3 / 2)` is real for some complex `f`
        return real(_LagrangianStrainToVolume(x.v0, f))
    else
        throw(DomainError("volume will be complex!"))
    end
end
_LagrangianStrainToVolume(v0, f) = v0 * (2f + 1)^_1½

struct NaturalStrainToVolume{T} <: StrainToVolume{T}
    v0::T
end
(x::NaturalStrainToVolume)(f) = _NaturalStrainToVolume(x.v0, f)
function (x::NaturalStrainToVolume)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(f)) / 60)  # `exp(3f)` is real for some complex `f`
        return real(_InfinitesimalStrainToVolume(x.v0, f))
    else
        throw(DomainError("volume will be complex!"))
    end
end
_NaturalStrainToVolume(v0, f) = v0 * exp(3f)

struct InfinitesimalStrainToVolume{T} <: StrainToVolume{T}
    v0::T
end
(x::InfinitesimalStrainToVolume)(f) = _InfinitesimalStrainToVolume(x.v0, f)
function (x::InfinitesimalStrainToVolume)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(1 - f)) / 60)  # `(1 - f)^3` is real for some complex `f`
        return real(_InfinitesimalStrainToVolume(x.v0, f))
    else
        throw(DomainError("volume will be complex!"))
    end
end
_InfinitesimalStrainToVolume(v0, f) = v0 / (1 - f)^3

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

end
