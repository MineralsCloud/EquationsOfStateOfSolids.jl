module EquationsOfStateOfSolids

using Unitful: AbstractQuantity, NoUnits, ħ, me

const FERMI_GAS_CONSTANT = (3π^2)^(2 / 3) * ħ^2 / 5 / me

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

Base.:^(x, ::TwoThirds) = x^(2 // 3)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::TwoThirds) = x^(2 / 3)
Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::TwoThirds) =
    isnegative(x) ? complex(x)^(2 / 3) : x^(2 / 3)
Base.:^(x, ::OneThird) = x^(1 // 3)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::OneThird) = x^(1 / 3)
Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::OneThird) =
    isnegative(x) ? complex(x)^(1 / 3) : x^(1 / 3)
Base.:^(x, ::FiveHalves) = x^(5 // 2)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::FiveHalves) = sqrt(x^5)
Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::FiveHalves) =
    isnegative(x) ? sqrt(complex(x^5)) : sqrt(x^5)
Base.:^(x, ::ThreeHalves) = x^(3 // 2)
Base.:^(x::Union{Complex,AbstractQuantity{<:Complex}}, ::ThreeHalves) = sqrt(x^3)
Base.:^(x::Union{Real,AbstractQuantity{<:Real}}, ::ThreeHalves) =
    isnegative(x) ? sqrt(complex(x^3)) : sqrt(x^3)

isnegative(x) = x < zero(x)  # Do not export!

include("FiniteStrains.jl")
include("collections/collections.jl")
include("vsolve.jl")
include("Fitting.jl")

end
