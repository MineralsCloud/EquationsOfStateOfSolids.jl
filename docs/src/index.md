```@meta
CurrentModule = EquationsOfStateOfSolids
```

# EquationsOfStateOfSolids

Documentation for [EquationsOfStateOfSolids](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code is [hosted on GitHub](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [@singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Package Features

This package implements some _equations of state_ (EOS) of solids which are
useful in research. It currently includes:

1. `Murnaghan1st` EOS
2. Birch–Murnaghan EOS family:
   1. `BirchMurnaghan2nd`
   2. `BirchMurnaghan3rd`
   3. `BirchMurnaghan4th`
3. `Vinet` EOS
4. Poirier–Tarantola EOS family:
   1. `PoirierTarantola2nd`
   2. `PoirierTarantola3rd`

The formulae are referenced from Ref. 1.

This package also includes linear and nonlinear fitting methods,
which are also referenced from Ref. 1.

- Calculate the energy, pressure, and bulk modulus of an `EquationOfStateOfSolid` on a
  volume (an array of volumes).
- Fit an `EquationOfStateOfSolid` on a series of ``E(V)`` data using the least-squares fitting
  method or a linear fitting method.
- Find the corresponding volume of energy, or pressure, given an `EquationOfStateOfSolid`.
- Handle unit conversion automatically in the above features.

The old [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl)
package has been superseded by `EquationsOfStateOfSolids.jl`.
So please just use `EquationsOfStateOfSolids.jl`.

### References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
2. [R. J. Angel, M. Alvaro, J. Gonzalez-Platas, _Zeitschrift Für Kristallographie - Cryst Mater_. **229**, 405–419 (2014).](https://www.degruyter.com/document/doi/10.1515/zkri-2013-1711/html)

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
pkg> add EquationsOfStateOfSolids
```

Or, equivalently, via the `Pkg` API:

```@repl
import Pkg; Pkg.add("EquationsOfStateOfSolids")
```

## Documentation

- [**STABLE**](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/stable) — **documentation of the most recently tagged version.**
- [**DEV**](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev) — _documentation of the in-development version._

## Project status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and contributions

Usage questions can be posted on
[our discussion page](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/discussions).

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues)
if you encounter any problems. The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual outline

```@contents
Pages = [
    "installation.md",
    "plotting.md",
    "interoperability.md",
    "portability.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
    "troubleshooting.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = [
    "api/collections.md",
    "api/finitestrains.md",
    "api/fitting.md",
]
```

### [Index](@id main-index)

```@index
Pages = [
    "api/collections.md",
    "api/finitestrains.md",
    "api/fitting.md",
]
```
