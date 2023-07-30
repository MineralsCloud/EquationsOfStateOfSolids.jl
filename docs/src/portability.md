# How to make your data portable?

After an equation-of-state-fitting, for instance, you want to save the data
to share with a colleague or for future use. Julia provides
several ways to do this. Below I will list one recommended way: saving it to
a JLD format by [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) package.

> JLD is a specific "dialect" of HDF5, a cross-platform, multi-language data storage format most frequently used for scientific data.

1. Install [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) and
   [FileIO.jl](https://github.com/JuliaIO/FileIO.jl) packages.

   ```@repl
   using Pkg
   Pkg.add("FileIO")
   Pkg.add("JLD2")
   ```

2. Create some `EquationsOfStateOfSolids`s:

   ```@repl
   using EquationsOfStateOfSolids, Unitful, UnitfulAtomic
   m = Murnaghan(224.501825, 0.00060479524074699499, 3.723835, -323.417686)
   bm = BirchMurnaghan3rd(224.4445656763778u"bohr^3", 9.194980249913018u"GPa", 3.7403684211716297, -161.70885710742223u"hartree")
   ```

3. Save them to file `"eos.jld2"`:

   ```@repl
   using JLD2, FileIO
   @save "eos.jld2" m bm
   ```

4. On another computer, or some days later, load them into REPL:

   ```@repl
   using EquationsOfStateOfSolids, Unitful, UnitfulAtomic
   @load "eos.jld2" m bm
   ```

   Now variables `m` and `bm` represent the original `Parameters`:

   ```@repl
   m.b0
   m.b′0
   bm.v0
   bm.b0
   ```

For more details on the JLD format, please refer to
[JLD.jl's doc](https://github.com/JuliaIO/JLD.jl/blob/master/doc/jld.md),
[JLD2.jl's doc](https://github.com/JuliaIO/JLD2.jl/blob/master/README.md) or
[this discussion](https://discourse.julialang.org/t/jld-jl-vs-jld2-jl/15287).
