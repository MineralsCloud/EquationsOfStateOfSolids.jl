module FiniteStrains

using ..EquationsOfStateOfSolids: _⅔, _⅓, _1½

abstract type FiniteStrain end
struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

"""
    ToEulerianStrain(v0)
    ToLagrangianStrain(v0)
    ToNaturalStrain(v0)
    ToInfinitesimalStrain(v0)

Calculate the finite strain of `v` based on the reference volume `v0`.

!!! info
    See the formulae on the [`Gibbs2` paper](https://www.sciencedirect.com/science/article/pii/S0010465511001470) Table 3.

# Examples
```jldoctest
julia> f = ToEulerianStrain(10);

julia> f(9)
0.036382991447572066

julia> f = ToEulerianStrain(100u"nm^3");

julia> f(90u"nm^3")
0.036382991447572066

julia> g = inv(f);

julia> g ∘ f == f ∘ g == identity
true
```
"""
struct To{S<:FiniteStrain,T}
    v0::T
end
To{S}(v0::T) where {S,T} = To{S,T}(v0)
const ToEulerianStrain = To{EulerianStrain}
const ToLagrangianStrain = To{LagrangianStrain}
const ToNaturalStrain = To{NaturalStrain}
const ToInfinitesimalStrain = To{InfinitesimalStrain}
(x::ToEulerianStrain)(v) = ((x.v0 / v)^_⅔ - 1) / 2
(x::ToLagrangianStrain)(v) = ((v / x.v0)^_⅔ - 1) / 2
(x::ToNaturalStrain)(v) = log(v / x.v0) / 3
(x::ToInfinitesimalStrain)(v) = 1 - (x.v0 / v)^_⅓

"""
    FromEulerianStrain(v0)
    FromLagrangianStrain(v0)
    FromNaturalStrain(v0)
    FromInfinitesimalStrain(v0)

Calculate the original volume `v` from the finite strain `f` based on the reference volume `v0`.

!!! info
    See the formulae on the [`Gibbs2` paper](https://www.sciencedirect.com/science/article/pii/S0010465511001470) Table 3.

# Examples
```jldoctest
julia> g = FromEulerianStrain(10);

julia> g(0.036382991447572066)
9.000000000000002

julia> g = FromEulerianStrain(100u"nm^3");

julia> g(0.036382991447572066)
90.00000000000001 nm³

julia> f = inv(g);

julia> f ∘ g == g ∘ f == identity
true
```
"""
struct From{S<:FiniteStrain,T}
    v0::T
end
From{S}(v0::T) where {S,T} = From{S,T}(v0)
const FromEulerianStrain = From{EulerianStrain}
const FromLagrangianStrain = From{LagrangianStrain}
const FromNaturalStrain = From{NaturalStrain}
const FromInfinitesimalStrain = From{InfinitesimalStrain}
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

Base.:∘(x::From{T}, y::To{T}) where {T} =
    x.v0 == y.v0 ? identity : error("undefined transformation!")
Base.:∘(x::To{T}, y::From{T}) where {T} = y ∘ x

Base.inv(x::From{T}) where {T} = To{T}(x.v0)
Base.inv(x::To{T}) where {T} = From{T}(x.v0)

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
