```@meta
CurrentModule = EquationsOfStateOfSolids
```

# EquationsOfStateOfSolids

## Package Features

- Calculate energy, pressure, and bulk modulus of an EOS `Parameters` on a (an)
  volume (array of volumes).
- Fit an `EquationOfState` on a series of volumes using least-squares fitting
  method.
- Fit an `EquationOfState` on a series of volumes linearly.
- Find the corresponding volume of an EOS `Parameters` given an (a) energy,
  pressure, and bulk modulus.
- Handle unit conversion automatically in the above features, take any unit.

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The old [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl)
package has been superseded by `EquationsOfStateOfSolids.jl`.
So please just use `EquationsOfStateOfSolids.jl`.

The code is [hosted on GitHub](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Compatibility

- [Julia version: `v1.0.0` to `v1.6.1`](https://julialang.org/downloads/)
- Dependencies:
  - [`AutoHashEquals.jl`](https://github.com/andrewcooke/AutoHashEquals.jl) `v0.2.0` and above
  - [`ConstructionBase.jl`](https://github.com/JuliaObjects/ConstructionBase.jl) `v1.0` and above
  - [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl) `v4.0.0` and above
  - [`LsqFit.jl`](https://github.com/JuliaNLSolvers/LsqFit.jl) `v0.8.0` and above
  - [`PolynomialRoots.jl`](https://github.com/giordano/PolynomialRoots.jl) `v1.0.0` and above
  - [`Polynomials.jl`](https://github.com/JuliaMath/Polynomials.jl) `v0.8.0` and above
  - [`Roots.jl`](https://github.com/JuliaMath/Roots.jl) `v0.8.0` and above
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
  - [`UnPack.jl`](https://github.com/mauro3/UnPack.jl) `v1.0.0` and above
  - [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) `v0.18.0` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Manual Outline

```@contents
Pages = [
    "installation.md",
    "develop.md",
    "portability.md",
    "interoperability.md",
    "plotting.md",
    "faq.md",
    "api/collections.md",
    "api/fitting.md",
]
Depth = 3
```

## [Index](@id main-index)

```@index
```
