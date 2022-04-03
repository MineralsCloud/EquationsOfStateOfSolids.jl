<div align="center">
  <img src="./docs/src/assets/logo.png" height="200"><br>
</div>

# EquationsOfStateOfSolids

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev)
[![Build Status](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/actions)
[![pipeline status](https://gitlab.com/singularitti/equationsofstateofsolids.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/equationsofstateofsolids.jl/-/pipelines)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/EquationsOfStateOfSolids.jl?svg=true)](https://ci.appveyor.com/project/singularitti/EquationsOfStateOfSolids-jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Coverage](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/E/EquationsOfStateOfSolids.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![GitHub license](https://img.shields.io/github/license/MineralsCloud/EquationsOfStateOfSolids.jl)](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/blob/master/LICENSE)

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

This package also includes linear and nonlinear fitting methods, also referenced
from Ref. 1.

See its
[documentation](https://mineralscloud.github.io/EquationsOfStateOfSolids.jl/stable).

## Compatibility

- [Julia version: above `v1.0.0`](https://julialang.org/downloads/)
- Dependencies:
  - [`StructHelpers.jl`](https://github.com/jw3126/StructHelpers.jl) `v0.1.0` and above
  - [`ConstructionBase.jl`](https://github.com/JuliaObjects/ConstructionBase.jl) `v1.0` and above
  - [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl) `v4.0.0` and above
  - [`LsqFit.jl`](https://github.com/JuliaNLSolvers/LsqFit.jl) `v0.8.0` and above
  - [`PolynomialRoots.jl`](https://github.com/giordano/PolynomialRoots.jl) `v1.0.0` and above
  - [`Polynomials.jl`](https://github.com/JuliaMath/Polynomials.jl) `v0.8.0` and above
  - [`Roots.jl`](https://github.com/JuliaMath/Roots.jl) `v0.8.0` and above
  - [`UnPack.jl`](https://github.com/mauro3/UnPack.jl) `v1.0.0` and above
  - [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) `v0.18.0` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
2. [R. J. Angel, M. Alvaro, J. Gonzalez-Platas, *Zeitschrift Für Kristallographie - Cryst Mater*. **229**, 405–419 (2014).](https://www.degruyter.com/document/doi/10.1515/zkri-2013-1711/html)
