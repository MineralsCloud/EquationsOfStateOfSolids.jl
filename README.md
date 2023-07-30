<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/EquationsOfStateOfSolids.jl/master/docs/src/assets/logo.png" height="200"><br>
</div>

# EquationsOfStateOfSolids

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                        **Others**                                         |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev
[gha-img]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/MineralsCloud/EquationsOfStateOfSolids.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/EquationsOfStateOfSolids-jl
[cirrus-img]: https://api.cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl
[gitlab-img]: https://gitlab.com/singularitti/EquationsOfStateOfSolids.jl/badges/main/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/EquationsOfStateOfSolids.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl
[license-img]: https://img.shields.io/github/license/MineralsCloud/EquationsOfStateOfSolids.jl
[license-url]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/blob/main/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

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
The formulae are referenced from [Ref. 1][Ref. 1].

- Calculate the energy, pressure, and bulk modulus of an `EquationOfStateOfSolid` on a
  volume (an array of volumes).
- Fit an `EquationOfStateOfSolid` on a series of ``E(V)`` data using the least-squares fitting
  method or a linear fitting method.
- Find the corresponding volume of energy, or pressure, given an `EquationOfStateOfSolid`.
- Handle unit conversion automatically in the above features.

The old [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl)
package has been superseded by `EquationsOfStateOfSolids.jl`.
So please just use `EquationsOfStateOfSolids.jl`.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add EquationsOfStateOfSolids
```

Or, equivalently, via the [`Pkg` API](https://pkgdocs.julialang.org/v1/getting-started/):

```julia
julia> import Pkg; Pkg.add("EquationsOfStateOfSolids")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

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

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=MineralsCloud/EquationsOfStateOfSolids.jl&type=Date)](https://star-history.com/#MineralsCloud/EquationsOfStateOfSolids.jl&Date)

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).][Ref. 1]
2. [R. J. Angel, M. Alvaro, J. Gonzalez-Platas, _Zeitschrift Für Kristallographie - Cryst Mater_. **229**, 405–419 (2014).](https://www.degruyter.com/document/doi/10.1515/zkri-2013-1711/html)

[Ref. 1]: https://www.sciencedirect.com/science/article/pii/S0010465511001470
