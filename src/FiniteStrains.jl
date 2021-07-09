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
    if f < -1 / 2
        @warn "strain $f < -0.5! Volume will be complex!"
    end
    return _FromEulerianStrain(x.v0, f)
end
function (x::FromEulerianStrain)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(2f + 1)) / 120)  # `(2f + 1)^(3 / 2)` is real for some complex `f`
        return real(_FromEulerianStrain(x.v0, f))
    else
        @warn "volume will be complex with strain $f."
        return _FromEulerianStrain(x.v0, f)
    end
end
_FromEulerianStrain(v0, f) = v0 / (2f + 1)^_1½
# Lagrangian strain
function (x::FromLagrangianStrain)(f)
    if f < -1 / 2
        @warn "strain $f < -0.5! Volume will be complex!"
    end
    return _FromLagrangianStrain(x.v0, f)
end
function (x::FromLagrangianStrain)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(2f + 1)) / 120)  # `(2f + 1)^(3 / 2)` is real for some complex `f`
        return real(_FromLagrangianStrain(x.v0, f))
    else
        @warn "volume will be complex with strain $f."
        return _FromLagrangianStrain(x.v0, f)
    end
end
_FromLagrangianStrain(v0, f) = v0 * (2f + 1)^_1½
# Natural strain
(x::FromNaturalStrain)(f) = x.v0 * exp(3f)
function (x::FromNaturalStrain)(f::Complex)
    v = x.v0 * exp(3f)
    if isreal(f) || isinteger(rad2deg(angle(f)) / 60)  # `exp(3f)` is real for some complex `f`
        return real(v)
    else
        @warn "volume will be complex with strain $f."
        return v
    end
end
(x::FromInfinitesimalStrain)(f) = _FromInfinitesimalStrain(x.v0, f)
function (x::FromInfinitesimalStrain)(f::Complex)
    if isreal(f)
        return x(real(f))
    elseif isinteger(rad2deg(angle(1 - f)) / 60)  # `(1 - f)^3` is real for some complex `f`
        return real(_FromInfinitesimalStrain(x.v0, f))
    else
        @warn "volume will be complex with strain $f."
        return _FromInfinitesimalStrain(x.v0, f)
    end
end
_FromInfinitesimalStrain(v0, f) = v0 / (1 - f)^3

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
