# Troubleshooting

```@contents
Pages = ["troubleshooting.md"]
Depth = 2
```

This page collects some possible errors you may encounter and trick how to fix them.
If you have some questions about how to use this code, you are welcome to
[discuss with us](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/discussions).

If you have additional tips, please either
[report an issue](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/issues/new) or
[submit a PR](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl/compare) with suggestions.

## Installation problems

### Cannot find the `julia` executable

Make sure you have Julia installed in your environment. Please download the latest
[stable version](https://julialang.org/downloads/#current_stable_release) for your platform.
If you are using a *nix system, the recommended way is to use
[Juliaup](https://github.com/JuliaLang/juliaup). If you do not want to install Juliaup
or you are using other platforms that Julia supports, download the corresponding binaries.
Then, create a symbolic link to the Julia executable. If the path is not in your `$PATH`
environment variable, export it to your `$PATH`.

Some clusters, like
[Habanero](https://confluence.columbia.edu/confluence/display/rcs/Habanero+HPC+Cluster+User+Documentation),
[Comet](https://www.sdsc.edu/support/user_guides/comet.html),
or [Expanse](https://www.sdsc.edu/services/hpc/expanse/index.html),
already have Julia installed as a module, you may
just `module load julia` to use it. If not, either install by yourself or contact your
administrator.

## Loading EquationsOfStateOfSolids

### Julia compiles/loads slow

First, we recommend you download the latest version of Julia. Usually, the newest version
has the best performance.

If you just want Julia to do a simple task and only once, you could start the Julia REPL with

```bash
julia --compile=min
```

to minimize compilation or

```bash
julia --optimize=0
```

to minimize optimizations, or just use both. Or you could make a system image
and run with

```bash
julia --sysimage custom-image.so
```

See [Fredrik Ekre's talk](https://youtu.be/IuwxE3m0_QQ?t=313) for details.

## How to make a `Vector` from a `Parameters`?

A suggested way is to use the
[`IterTools.fieldvalues` function](https://juliacollections.github.io/IterTools.jl/latest/index.html#IterTools.fieldvalues):

```@repl
using IterTools
eos = BirchMurnaghan4th(1, 2.0, 3, 4)
collect(fieldvalues(eos))
```

It is lazy and fast.

Or, write a non-lazy version of `fieldvalues` manually:

```@repl
fieldvalues(eos::EquationOfState) = [getfield(eos, i) for i in 1:nfields(eos)]
fieldvalues(eos)
```

## `linfit` does not work with `BigFloat`?

`LinearAlgebra` does not support SVD for matrices with `BigFloat`
elements by default. You need to install
[`GenericSVD.jl`](https://github.com/JuliaLinearAlgebra/GenericSVD.jl) first
and then `using GenericSVD`. Then it should work.
