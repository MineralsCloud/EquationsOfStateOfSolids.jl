<div align="center">
  <img src="./docs/src/assets/logo.png" height="200"><br>
</div>

# EquationsOfStateOfSolids

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev)
[![Build Status](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/EquationsOfStateOfSolids.jl?svg=true)](https://ci.appveyor.com/project/singularitti/EquationsOfStateOfSolids-jl)
[![Build Status](https://cloud.drone.io/api/badges/MineralsCloud/EquationsOfStateOfSolids.jl/status.svg)](https://cloud.drone.io/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Coverage](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/EquationsOfStateOfSolids.jl)
[![Coverage](https://coveralls.io/repos/github/MineralsCloud/EquationsOfStateOfSolids.jl/badge.svg?branch=master)](https://coveralls.io/github/MineralsCloud/EquationsOfStateOfSolids.jl?branch=master)

This package implements some _equations of state_ (EOS) of solids which are
useful in research. It currently includes:

1. `Murnaghan` EOS
2. Birch–Murnaghan EOS family:
   1. `BirchMurnaghan2nd`
   2. `BirchMurnaghan3rd`
   3. `BirchMurnaghan4th`
3. `Vinet` EOS
4. Poirier–Tarantola EOS family:
   1. `PoirierTarantola2nd`
   2. `PoirierTarantola3rd`
   3. `PoirierTarantola4th`
5. `AntonSchmidt` EOS (experimental)

The formulae are referenced from Ref. 1.

This package also includes linear and nonlinear fitting methods, also referenced
from Ref. 1.

See its
[documentation](https://mineralscloud.github.io/EquationsOfStateOfSolids.jl/stable).

## Compatibility

- [Julia version: `v1.0.0` and above](https://julialang.org/downloads/)
- Dependencies: see `Project.toml` [`deps` field](Project.toml#L7-L18) and
  [`compat` field](Project.toml#L21-L31)
- OS: macOS, Linux, and Windows

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
