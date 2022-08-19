<div align="center">
  <img src="https://raw.githubusercontent.com/MineralsCloud/EquationsOfStateOfSolids.jl/master/docs/src/assets/logo.png" height="200"><br>
</div>

# EquationsOfStateOfSolids

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                        **Others**                                         |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://singularitti.github.io/MyPkgTemplates.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://singularitti.github.io/MyPkgTemplates.jl/dev
[gha-img]: https://github.com/singularitti/MyPkgTemplates.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/singularitti/MyPkgTemplates.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/singularitti/MyPkgTemplates.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/MyPkgTemplates-jl
[cirrus-img]: https://api.cirrus-ci.com/github/singularitti/MyPkgTemplates.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/singularitti/MyPkgTemplates.jl
[gitlab-img]: https://gitlab.com/singularitti/MyPkgTemplates.jl/badges/master/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/MyPkgTemplates.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/singularitti/MyPkgTemplates.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/singularitti/MyPkgTemplates.jl
[license-img]: https://img.shields.io/github/license/singularitti/MyPkgTemplates.jl
[license-url]: https://github.com/singularitti/MyPkgTemplates.jl/blob/master/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle

The code is [hosted on GitHub](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [@singularitti](https://github.com/singularitti).
You are very welcome to contribute.

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

The formulae are referenced from Ref. 1.

This package also includes linear and nonlinear fitting methods,
which are also referenced from Ref. 1.

- Calculate the energy, pressure, and bulk modulus of an `EquationOfStateOfSolid` on a
  volume (an array of volumes).
- Fit an `EquationOfStateOfSolid` on a series of `E(V)` data using the least-squares fitting
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

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and contributions

Usage questions can be posted on [our discussion page][discussions-url].

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems. The [contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

[discussions-url]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/discussions
[issues-url]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues
[contrib-url]: https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/discussions

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
2. [R. J. Angel, M. Alvaro, J. Gonzalez-Platas, _Zeitschrift Für Kristallographie - Cryst Mater_. **229**, 405–419 (2014).](https://www.degruyter.com/document/doi/10.1515/zkri-2013-1711/html)
