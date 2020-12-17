# Plotting

Package
[`EquationOfStateRecipes.jl`](https://github.com/MineralsCloud/EquationOfStateRecipes.jl)
provides some default themes for plotting an `EquationOfStateOfSolids`. First,
try to install the [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) package
by

```julia
julia> using Pkg

julia> Pkg.add("Plots")

julia> using Plots
```

Then install `EquationOfStateRecipes.jl` with

```julia
julia> pkg"add https://github.com/MineralsCloud/EquationOfStateRecipes.jl.git"

julia> using EquationOfStateRecipes
```

Finally, load `EquationsOfStateOfSolids.jl` and build your own
`EquationOfStateOfSolids`:

```julia
julia> using EquationsOfStateOfSolids.Collections

julia> eos = EnergyEos(Murnaghan(224.501825u"bohr^3", 8.896845u"GPa", 3.723835, -323.417686u"Ry"));

julia> plot(eos)

julia> plot!(eos, (0.8:0.01:1.2) * eos.param.v0)

julia> scatter!(eos, (0.5:0.1:1) * eos.param.v0)
```

![fig](https://i.loli.net/2020/12/16/BrbLsZlmKvy6hTi.png)

Have fun!
