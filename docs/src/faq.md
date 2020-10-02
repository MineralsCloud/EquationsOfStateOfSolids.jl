# FAQ

## How to make a `Vector` from a `Parameters`?

A suggested way is to use the
[`IterTools.fieldvalues` function](https://juliacollections.github.io/IterTools.jl/latest/index.html#IterTools.fieldvalues):

```julia
julia> using IterTools

julia> eos = BirchMurnaghan4th(1, 2.0, 3, 4)
BirchMurnaghan4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> collect(fieldvalues(eos))
5-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 0.0
```

It is lazy and fast.

Or to write a non-lazy version of `fieldvalues` manually:

```julia
julia> fieldvalues(eos::EquationOfState) = [getfield(eos, i) for i in 1:nfields(eos)]
fieldvalues (generic function with 1 method)

julia> fieldvalues(eos)
5-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 0.0
```

It is slower than `IterTools.fieldvalues`. Use it with care.

## `linfit` does not work with `BigFloat`?

`LinearAlgebra` by default does not support SVD for matrices with `BigFloat` elements.
You need to install [`GenericSVD.jl`](https://github.com/JuliaLinearAlgebra/GenericSVD.jl)
first then `using GenericSVD`. And then it should work.
