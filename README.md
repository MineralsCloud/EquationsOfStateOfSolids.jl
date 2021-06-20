<div align="center">
  <img src="./docs/src/assets/logo.png" height="200"><br>
</div>

# EquationsOfStateOfSolids

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev)
[![Build Status](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/actions)
[![pipeline status](https://gitlab.com/singularitti/equationsofstateofsolids.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/equationsofstateofsolids.jl/-/pipelines)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/EquationsOfStateOfSolids.jl?svg=true)](https://ci.appveyor.com/project/singularitti/EquationsOfStateOfSolids-jl)
[![Build Status](https://cloud.drone.io/api/badges/MineralsCloud/EquationsOfStateOfSolids.jl/status.svg)](https://cloud.drone.io/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Coverage](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Coverage](https://coveralls.io/repos/github/MineralsCloud/EquationsOfStateOfSolids.jl/badge.svg?branch=master)](https://coveralls.io/github/MineralsCloud/EquationsOfStateOfSolids.jl?branch=master)

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

- [Julia version: `v1.0.0` to `v1.6.1`](https://julialang.org/downloads/)
- Dependencies:
  - [`AutoHashEquals.jl`](https://github.com/andrewcooke/AutoHashEquals.jl) `v0.2.0` and above
  - [`Compat.jl`](https://github.com/JuliaLang/Compat.jl) `v3.1.0` and above
  - [`Configurations.jl`](https://github.com/Roger-luo/Configurations.jl) `v0.3.0` and above
  - [`ConstructionBase.jl`](https://github.com/JuliaObjects/ConstructionBase.jl) `v1.0` and above
  - [`EquationsOfState.jl`](https://github.com/MineralsCloud/EquationsOfState.jl) `v4.0.0` and above
  - [`LsqFit.jl`](https://github.com/JuliaNLSolvers/LsqFit.jl) `v0.8.0` and above
  - [`PolynomialRoots.jl`](https://github.com/giordano/PolynomialRoots.jl) `v1.0.0` and above
  - [`Polynomials.jl`](https://github.com/JuliaMath/Polynomials.jl) `v0.8.0` and above
  - [`Roots.jl`](https://github.com/JuliaMath/Roots.jl) `v0.8.0` and above
  - [`Roots.jl`](https://github.com/JuliaMath/Roots.jl) `v0.8.0` and above
  - [`UnPack.jl`](https://github.com/mauro3/UnPack.jl) `v1.0.0` and above
  - [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) `v0.18.0` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
