```@meta
CurrentModule = EquationsOfStateOfSolids
```

# EquationsOfStateOfSolids

Documentation for [EquationsOfStateOfSolids](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

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

This package also includes linear and nonlinear fitting methods.
The formulae are referenced from [Ref. 1](https://www.sciencedirect.com/science/article/pii/S0010465511001470).

- Calculate the energy, pressure, and bulk modulus of an `EquationOfStateOfSolid` on a
  volume (an array of volumes).
- Fit an `EquationOfStateOfSolid` on a series of ``E(V)`` data using the least-squares fitting
  method or a linear fitting method.
- Find the corresponding volume of energy, or pressure, given an `EquationOfStateOfSolid`.
- Handle unit conversion automatically in the above features.

The old [EquationsOfState.jl](https://github.com/MineralsCloud/EquationsOfState.jl)
package has been superseded by EquationsOfStateOfSolids.jl.
So please just use EquationsOfStateOfSolids.jl.

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

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual outline

```@contents
Pages = [
    "man/installation.md",
    "man/definitions.md",
    "man/examples.md",
    "man/troubleshooting.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
