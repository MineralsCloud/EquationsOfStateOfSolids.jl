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
struct StrainFromVolume{S<:FiniteStrain} end
const EulerianStrainFromVolume = StrainFromVolume{Eulerian}
const LagrangianStrainFromVolume = StrainFromVolume{Lagrangian}
const NaturalStrainFromVolume = StrainFromVolume{Natural}
const InfinitesimalStrainFromVolume = StrainFromVolume{Infinitesimal}
struct StrainFromVolumeWithReferenceVolume{S<:FiniteStrain,T<:Number}
    v0::T
end
function StrainFromVolumeWithReferenceVolume{S}(v0) where {S}
    return StrainFromVolumeWithReferenceVolume{S,typeof(v0)}(v0)
end
StrainFromVolume(::T) where {T} = StrainFromVolume{T}()
(::StrainFromVolume{T})(v0) where {T} = StrainFromVolumeWithReferenceVolume{T}(v0)
(x::StrainFromVolumeWithReferenceVolume{Eulerian})(v) = ((x.v0 / v)^_⅔ - 1) / 2
(x::StrainFromVolumeWithReferenceVolume{Lagrangian})(v) = ((v / x.v0)^_⅔ - 1) / 2
(x::StrainFromVolumeWithReferenceVolume{Natural})(v) = log(v / x.v0) / 3
(x::StrainFromVolumeWithReferenceVolume{Infinitesimal})(v) = 1 - (x.v0 / v)^_⅓

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
struct VolumeFromStrain{S<:FiniteStrain} end
const VolumeFromEulerianStrain = VolumeFromStrain{Eulerian}
const VolumeFromLagrangianStrain = VolumeFromStrain{Lagrangian}
const VolumeFromNaturalStrain = VolumeFromStrain{Natural}
const VolumeFromInfinitesimalStrain = VolumeFromStrain{Infinitesimal}
struct VolumeFromStrainWithReferenceVolume{S<:FiniteStrain,T<:Number}
    v0::T
end
function VolumeFromStrainWithReferenceVolume{S}(v0) where {S}
    return VolumeFromStrainWithReferenceVolume{S,typeof(v0)}(v0)
end
VolumeFromStrain(::T) where {T} = VolumeFromStrain{T}()
(::VolumeFromStrain{T})(v0) where {T} = VolumeFromStrainWithReferenceVolume{T}(v0)
function (x::VolumeFromStrainWithReferenceVolume{Eulerian})(f)
    v = x.v0 / (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
function (x::VolumeFromStrainWithReferenceVolume{Lagrangian})(f)
    v = x.v0 * (2f + 1)^_1½
    return isreal(v) ? real(v) : v
end
function (x::VolumeFromStrainWithReferenceVolume{Natural})(f)
    v = x.v0 * exp(3f)
    return isreal(v) ? real(v) : v
end
function (x::VolumeFromStrainWithReferenceVolume{Infinitesimal})(f)
    v = x.v0 / (1 - f)^3
    return isreal(v) ? real(v) : v
end

function Base.:∘(
    x::VolumeFromStrainWithReferenceVolume{T}, y::StrainFromVolumeWithReferenceVolume{T}
) where {T}
    return x.v0 == y.v0 ? identity : error("undefined transformation!")
end
function Base.:∘(
    x::StrainFromVolumeWithReferenceVolume{T}, y::VolumeFromStrainWithReferenceVolume{T}
) where {T}
    return y ∘ x
end

function Base.inv(x::VolumeFromStrainWithReferenceVolume{T}) where {T}
    return StrainFromVolumeWithReferenceVolume{T}(x.v0)
end
function Base.inv(x::StrainFromVolumeWithReferenceVolume{T}) where {T}
    return VolumeFromStrainWithReferenceVolume{T}(x.v0)
end

"""
    DerivativeOfStrain{Eulerian}(N)
    DerivativeOfStrain{Eulerian}(N)
    DerivativeOfStrain{Eulerian}(N)
    DerivativeOfStrain{Eulerian}(N)

Return a function of `v₀` and `v` that calculates the `N`th derivative of the strain
with respect to volume at `v` with reference volume `v₀`.
"""
struct DerivativeOfStrain{S<:FiniteStrain,N} end
DerivativeOfStrain{S}(N::Integer) where {S} = DerivativeOfStrain{S,UInt8(N)}()
(derivative::DerivativeOfStrain{Eulerian})(v₀) = v -> derivative(v₀, v)
(::DerivativeOfStrain{Eulerian,1})(v₀, v) = -(v₀ / v)^_⅔ / 3v
function (::DerivativeOfStrain{Eulerian,N})(v₀, v) where {N}
    return (1 - 3N) / 3v * DerivativeOfStrain{Eulerian}(N - 1)(v₀, v)
end
(::DerivativeOfStrain{Lagrangian,1})(v₀, v) = -(v / v₀)^_⅔ / 3v
function (::DerivativeOfStrain{Lagrangian,N})(v₀, v) where {N}
    return (5 - 3N) / 3v * DerivativeOfStrain{Lagrangian}(N - 1)(v₀, v)
end
(::DerivativeOfStrain{Natural,1})(::Number, v) = 1 / 3v
function (::DerivativeOfStrain{Natural,N})(v₀, v) where {N}
    return (1 - N) / v * DerivativeOfStrain{Natural}(N - 1)(v₀, v)
end
function (::DerivativeOfStrain{Infinitesimal,1})(v₀, v)
    return (1 - InfinitesimalStrainFromVolume(v₀)(v))^4 / 3 / v₀
end
function (::DerivativeOfStrain{Infinitesimal,N})(v₀, v) where {N}
    return (2 - 3N) / 3v * DerivativeOfStrain{Infinitesimal}(N - 1)(v₀, v)
end

end
