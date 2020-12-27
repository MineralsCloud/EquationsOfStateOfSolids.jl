module FiniteStrains

using Unitful: Volume

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
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end
const Eulerian = EulerianStrain
const Lagrangian = LagrangianStrain
const Natural = NaturalStrain
const Infinitesimal = InfinitesimalStrain
const V0 = Union{Real,Volume{<:Real}}

"""
    volume2strain(::EulerianStrain, v0)
    volume2strain(::LagrangianStrain, v0)
    volume2strain(::NaturalStrain, v0)
    volume2strain(::InfinitesimalStrain, v0)

Return a function of `v` that calculates the `FiniteStrain` from `v0`.

!!! info
    See the formulae on Ref. 1 Table 3.
"""
volume2strain(::Eulerian, v0::V0) = v -> ((v0 / v)^_⅔ - 1) / 2
volume2strain(::Lagrangian, v0::V0) = v -> ((v / v0)^_⅔ - 1) / 2
volume2strain(::Natural, v0::V0) = v -> log(v / v0) / 3
volume2strain(::Infinitesimal, v0::V0) = v -> 1 - (v0 / v)^_⅓

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

end
