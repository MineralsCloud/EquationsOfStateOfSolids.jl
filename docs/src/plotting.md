# Plotting

Package
[`EquationOfStateRecipes.jl`](https://github.com/MineralsCloud/EquationOfStateRecipes.jl)
provides some default themes for plotting an `EquationOfStateOfSolids`.
First, try to install the [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) package by

```@repl
using Pkg
Pkg.add("Plots")
```

Then install `EquationOfStateRecipes.jl` with

```@repl
Pkg.add("EquationOfStateRecipes")
```

Finally, load `EquationsOfStateOfSolids.jl` and plot:

```@repl
using EquationsOfStateOfSolids, Plots
eos = EnergyEquation(Murnaghan(224.501825u"bohr^3", 8.896845u"GPa", 3.723835, -323.417686u"Ry"))
plot(eos)
plot!(eos, (0.8:0.01:1.2) * eos.param.v0)
scatter!(eos, (0.5:0.1:1) * eos.param.v0)
```

![fig](https://i.loli.net/2020/12/16/BrbLsZlmKvy6hTi.png)

Have fun!
