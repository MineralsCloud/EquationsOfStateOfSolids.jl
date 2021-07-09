module FiniteStrains

using ..EquationsOfStateOfSolids: _⅔, _⅓, _1½

abstract type FiniteStrain end
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

struct ToStrain{S<:FiniteStrain,T}
    v0::T
end
ToStrain{S}(v0::T) where {S,T} = ToStrain{S,T}(v0)
const ToEulerianStrain = ToStrain{EulerianStrain}
const ToLagrangianStrain = ToStrain{LagrangianStrain}
const ToNaturalStrain = ToStrain{NaturalStrain}
const ToInfinitesimalStrain = ToStrain{InfinitesimalStrain}
(x::ToEulerianStrain)(v) = ((x.v0 / v)^_⅔ - 1) / 2
(x::ToLagrangianStrain)(v) = ((v / x.v0)^_⅔ - 1) / 2
(x::ToNaturalStrain)(v) = log(v / x.v0) / 3
(x::ToInfinitesimalStrain)(v) = 1 - (x.v0 / v)^_⅓

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
struct FromStrain{S<:FiniteStrain,T}
    v0::T
end
FromStrain{S}(v0::T) where {S,T} = FromStrain{S,T}(v0)
const FromEulerianStrain = FromStrain{EulerianStrain}
const FromLagrangianStrain = FromStrain{LagrangianStrain}
const FromNaturalStrain = FromStrain{NaturalStrain}
const FromInfinitesimalStrain = FromStrain{InfinitesimalStrain}
# Eulerian strain
function (x::FromEulerianStrain)(f)
    v = x.v0 / (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
# Lagrangian strain
function (x::FromLagrangianStrain)(f)
    v = x.v0 * (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
# Natural strain
function (x::FromNaturalStrain)(f)
    v = x.v0 * exp(3f)
    return isreal(v) ? real(v) : v
end
# Infinitesimal strain
function (x::FromInfinitesimalStrain)(f)
    v = x.v0 / (1 - f)^3
    return isreal(v) ? real(v) : v
end

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
            return (1 - ToInfinitesimalStrain(v0)(v))^4 / 3 / v0
        else  # Recursion
            return -(3deg - 2) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end

function straintype end

end
