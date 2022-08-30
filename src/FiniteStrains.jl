module FiniteStrains

using ..EquationsOfStateOfSolids: _⅔, _⅓, _1½

export EulerianStrainFromVolume,
    LagrangianStrainFromVolume,
    NaturalStrainFromVolume,
    InfinitesimalStrainFromVolume,
    VolumeFromEulerianStrain,
    VolumeFromLagrangianStrain,
    VolumeFromNaturalStrain,
    VolumeFromInfinitesimalStrain

abstract type FiniteStrain end
struct Eulerian <: FiniteStrain end
struct Lagrangian <: FiniteStrain end
struct Natural <: FiniteStrain end
struct Infinitesimal <: FiniteStrain end

"""
    EulerianStrainFromVolume(v0)
    LagrangianStrainFromVolume(v0)
    NaturalStrainFromVolume(v0)
    InfinitesimalStrainFromVolume(v0)

Calculate the finite strain of `v` based on the reference volume `v0`.

# Examples
```jldoctest
julia> f = EulerianStrainFromVolume(10);

julia> f(9)
0.036382991447572066

julia> f = EulerianStrainFromVolume(100u"nm^3");

julia> f(90u"nm^3")
0.036382991447572066

julia> g = inv(f);

julia> g ∘ f == f ∘ g == identity
true
```
"""
struct StrainFromVolume{S<:FiniteStrain,T<:Number}
    v0::T
end
StrainFromVolume{S}(v0::T) where {S,T} = StrainFromVolume{S,T}(v0)
const EulerianStrainFromVolume = StrainFromVolume{Eulerian}
const LagrangianStrainFromVolume = StrainFromVolume{Lagrangian}
const NaturalStrainFromVolume = StrainFromVolume{Natural}
const InfinitesimalStrainFromVolume = StrainFromVolume{Infinitesimal}
(x::EulerianStrainFromVolume)(v) = ((x.v0 / v)^_⅔ - 1) / 2
(x::LagrangianStrainFromVolume)(v) = ((v / x.v0)^_⅔ - 1) / 2
(x::NaturalStrainFromVolume)(v) = log(v / x.v0) / 3
(x::InfinitesimalStrainFromVolume)(v) = 1 - (x.v0 / v)^_⅓

"""
    VolumeFromEulerianStrain(v0)
    VolumeFromLagrangianStrain(v0)
    VolumeFromNaturalStrain(v0)
    VolumeFromInfinitesimalStrain(v0)

Calculate the original volume `v` from the finite strain `f` based on the reference volume `v0`.

# Examples
```jldoctest
julia> g = VolumeFromEulerianStrain(10);

julia> g(0.036382991447572066)
9.000000000000002

julia> g = VolumeFromEulerianStrain(100u"nm^3");

julia> g(0.036382991447572066)
90.00000000000001 nm³

julia> f = inv(g);

julia> f ∘ g == g ∘ f == identity
true
```
"""
struct VolumeFromStrain{S<:FiniteStrain,T<:Number}
    v0::T
end
VolumeFromStrain{S}(v0::T) where {S,T} = VolumeFromStrain{S,T}(v0)
const VolumeFromEulerianStrain = VolumeFromStrain{Eulerian}
const VolumeFromLagrangianStrain = VolumeFromStrain{Lagrangian}
const VolumeFromNaturalStrain = VolumeFromStrain{Natural}
const VolumeFromInfinitesimalStrain = VolumeFromStrain{Infinitesimal}
# Eulerian strain
function (x::VolumeFromEulerianStrain)(f)
    v = x.v0 / (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
# Lagrangian strain
function (x::VolumeFromLagrangianStrain)(f)
    v = x.v0 * (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
# Natural strain
function (x::VolumeFromNaturalStrain)(f)
    v = x.v0 * exp(3f)
    return isreal(v) ? real(v) : v
end
# Infinitesimal strain
function (x::VolumeFromInfinitesimalStrain)(f)
    v = x.v0 / (1 - f)^3
    return isreal(v) ? real(v) : v
end

function Base.:∘(x::VolumeFromStrain{T}, y::StrainFromVolume{T}) where {T}
    return x.v0 == y.v0 ? identity : error("undefined transformation!")
end
Base.:∘(x::StrainFromVolume{T}, y::VolumeFromStrain{T}) where {T} = y ∘ x

Base.inv(x::VolumeFromStrain{T}) where {T} = StrainFromVolume{T}(x.v0)
Base.inv(x::StrainFromVolume{T}) where {T} = VolumeFromStrain{T}(x.v0)

"""
    Dⁿᵥf(s::EulerianStrain, deg, v0)
    Dⁿᵥf(s::LagrangianStrain, deg, v0)
    Dⁿᵥf(s::NaturalStrain, deg, v0)
    Dⁿᵥf(s::InfinitesimalStrain, deg, v0)

Return a function of `v` that calculates the `deg`th order derivative of strain wrt volume from `v0`.
"""
function Dⁿᵥf(s::Eulerian, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v0 / v)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 1) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::Lagrangian, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return -(v / v0)^_⅔ / 3 / v
        else  # Recursion
            return -(3deg - 5) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::Natural, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return 1 / 3 / v
        else  # Recursion
            return -(deg - 1) / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end
function Dⁿᵥf(s::Infinitesimal, deg, v0)
    function (v)
        if isone(deg)  # Stop recursion
            return (1 - InfinitesimalStrainFromVolume(v0)(v))^4 / 3 / v0
        else  # Recursion
            return -(3deg - 2) / 3 / v * Dⁿᵥf(s, deg - 1, v0)(v)
        end
    end
end

end
