module EquationsOfStateOfSolids

using Unitful: AbstractQuantity, NoUnits

# See https://discourse.julialang.org/t/is-there-a-way-to-include-in-function-name/45378/3
abstract type Power end
struct TwoThirds <: Power end
struct OneThird <: Power end
struct FiveHalves <: Power end
struct ThreeHalves <: Power end

const _⅔ = TwoThirds()
const _⅓ = OneThird()
const _2½ = FiveHalves()
const _1½ = ThreeHalves()

Base.:^(x, ::TwoThirds) = x^(2//3)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::TwoThirds) = x^(2 / 3)
function Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::TwoThirds)
    return isnegative(x) ? complex(x)^(2 / 3) : x^(2 / 3)
end
Base.:^(x, ::OneThird) = x^(1//3)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::OneThird) = x^(1 / 3)
function Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::OneThird)
    return isnegative(x) ? complex(x)^(1 / 3) : x^(1 / 3)
end
Base.:^(x, ::FiveHalves) = x^(5//2)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::FiveHalves) = sqrt(x^5)
function Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::FiveHalves)
    return isnegative(x) ? sqrt(complex(x^5)) : sqrt(x^5)
end
Base.:^(x, ::ThreeHalves) = x^(3//2)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::ThreeHalves) = sqrt(x^3)
function Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::ThreeHalves)
    return isnegative(x) ? sqrt(complex(x^3)) : sqrt(x^3)
end

isnegative(x) = x < zero(x)  # Do not export!

include("FiniteStrains.jl")
include("collections/collections.jl")
include("vsolve.jl")
include("Fitting/Fitting.jl")

end
