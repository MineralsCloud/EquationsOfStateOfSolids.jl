# How to contribute

## Download the project

Similar to section "[Installation](@ref)", run

```@repl
using Pkg
pkg"dev EquationsOfStateOfSolids"
```

in Julia REPL.

Then the package will be cloned to your local machine at a path. On macOS, by default is
located at `~/.julia/dev/EquationsOfStateOfSolids` unless you modify the `JULIA_DEPOT_PATH`
environment variable. (See [Julia's official documentation](http://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_DEPOT_PATH-1)
on how to do this.) In the following text, we will call it `PKGROOT`.

## [Instantiate the project](@id instantiating)

Go to `PKGROOT`, start a new Julia session and run

```@repl
using Pkg
Pkg.instantiate()
```

## How to build docs

Usually, the up-to-state doc is available in
[here](https://MineralsCloud.github.io/EquationsOfStateOfSolids.jl/dev), but there are cases
where users need to build the doc themselves.

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
julia --color=yes docs/make.jl
```

in your terminal. In a while a folder `PKGROOT/docs/build` will appear. Open
`PKGROOT/docs/build/index.html` with your favorite browser and have fun!

## How to run tests

After [instantiating](@ref) the project, go to `PKGROOT`, run (without the `$` prompt)

```bash
julia --color=yes test/runtests.jl
```

in your terminal.
