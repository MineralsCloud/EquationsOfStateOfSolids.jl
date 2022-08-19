```@meta
CurrentModule = EquationsOfStateOfSolids
DocTestSetup  = :(using EquationsOfStateOfSolids)
```

# Collections

```@contents
Pages = ["collections.md"]
Depth = 3
```

The current `EquationOfStateOfSolidsParameters` are

```@repl
using TypeTree
tt(EquationsOfStateOfSolids.EquationOfStateOfSolidsParameters)
```

Here the leaves of the type tree are concrete types and can be constructed.

## Usage

### Construct a `EquationOfStateOfSolidsParameters` instance

We will use `BirchMurnaghan3rd` as an example.

A `BirchMurnaghan3rd` can be constructed from scratch, as shown above. It can
also be constructed from an existing `BirchMurnaghan3rd`, with
[`Setfield.jl`](https://github.com/jw3126/Setfield.jl)
[`@set!`](https://jw3126.github.io/Setfield.jl/stable/#Setfield.@set!-Tuple{Any})
macro:

```@repl
using Setfield
eos = Murnaghan1st(1, 2, 3.0)
@set! eos.v0 = 4
eos
```

To modify multiple fields (say, `:v0`, `:b′0`, `:b″0`, `:e0`) at a time, use
[`@batchlens`](https://tkf.github.io/Kaleido.jl/stable/#Kaleido.@batchlens) from
[`Kaleido.jl`](https://github.com/tkf/Kaleido.jl):

```@repl
using Setfield, Kaleido
lens = @batchlens(begin
           _.v0
           _.b′0
           _.b″0
           _.e0
       end)
eos = BirchMurnaghan4th(1, 2.0, 3, 4)
set(eos, lens, (5, 6, 7, 8))
```

Users can access `BirchMurnaghan3rd`'s elements by "dot notation":

```@repl
eos = BirchMurnaghan3rd(1, 2, 3, 4.0)
eos.v0
```

### Evaluate energy

The ``E(V)`` relation of equations of state are listed as below:

1. `Murnaghan`:

2. `BirchMurnaghan2nd`:

   ```math
   ```

3. `BirchMurnaghan3rd`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

   where ``x = V / V_0``, and
   ``f = \frac{ 1 }{ 2 } \bigg[ \bigg( \frac{ V_0 }{ V } \bigg)^{2/3} - 1 \bigg]``.

4. `BirchMurnaghan4th`:

   ```math
   E(V) = E_{0}+\frac{3}{8} V_{0} B_{0} f^{2}\left[\left(9 H-63 B_{0}^{\prime}+143\right) f^{2}+12\left(B_{0}^{\prime}-4\right) f+12\right].
   ```

   where ``H = B_0 B_0'' + (B_0')^2``.

5. `PoirierTarantola2nd`:

   ```math
   E(V) = E_{0}+\frac{1}{2} B_{0} V_{0} \ln ^{2} x.
   ```

6. `PoirierTarantola3rd`:

   ```math
   E(V) = E_{0}+\frac{1}{6} B_{0} V_{0} \ln ^{2} x\left[\left(B_{0}^{\prime}+2\right) \ln x+3\right].
   ```

7. `PoirierTarantola4th`:

   ```math
   E(V) = E_{0}+\frac{1}{24} B_{0} V_{0} \ln ^{2} x\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{2} x\right. \left.+4\left(B_{0}^{\prime}+2\right) \ln x+12\right\}.
   ```

   where ``H = B_0 B_0'' + (B_0')^2``.

8. `Vinet`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

9. `AntonSchmidt`:

   ```math
   E(V)=\frac{\beta V_{0}}{n+1}\left(\frac{V}{V_{0}}\right)^{n+1}\left[\ln \left(\frac{V}{V_{0}}\right)-\frac{1}{n+1}\right]+E_{\infty}.
   ```

### Evaluate pressure

The $P(V)$ relation of equations of state are listed as below:

1. `Murnaghan`:

   ```math
   1
   ```

2. `BirchMurnaghan2nd`:

   ```math
   P(V) = \frac{3}{2} B_{0}\left(x^{-7 / 3}-x^{-5 / 3}\right).
   ```

3. `BirchMurnaghan3rd`:

   ```math
   P(V) = \frac{3}{8} B_{0} \frac{x^{2 / 3}-1}{x^{10 / 3}}\left\{3 B_{0}^{\prime} x-16 x-3 x^{1 / 3}\left(B_{0}^{\prime}-4\right)\right\}.
   ```

4. `BirchMurnaghan4th`:

   ```math
   P(V) = \frac{1}{2} B_{0}(2 f+1)^{5 / 2}\left\{\left(9 H-63 B_{0}^{\prime}+143\right) f^{2}\right.\left.+9\left(B_{0}^{\prime}-4\right) f+6\right\}.
   ```

5. `PoirierTarantola2nd`:

   ```math
   P(V) = -\frac{B_{0}}{x} \ln x.
   ```

6. `PoirierTarantola3rd`:

   ```math
   P(V) = -\frac{B_{0} \ln x}{2 x}\left[\left(B_{0}^{\prime}+2\right) \ln x+2\right].
   ```

7. `PoirierTarantola4th`:

   ```math
   P(V) = -\frac{B_{0} \ln x}{6 x}\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{2} x+3\left(B_{0}^{\prime}+6\right) \ln x+6\right\}.
   ```

8. `Vinet`:

   ```math
   P(V) = 3 B_{0} \frac{1-\eta}{\eta^{2}} \exp \left\{-\frac{3}{2}\left(B_{0}^{\prime}-1\right)(\eta-1)\right\}.
   ```

9. `AntonSchmidt`:

   ```math
   P(V) = -\beta\left(\frac{V}{V_{0}}\right)^{n} \ln \left(\frac{V}{V_{0}}\right).
   ```

### Evaluate bulk modulus

The $B(V)$ relation of equations of state are listed as below:

1. `BirchMurnaghan2nd`:

   ```math
   B(V) = B_{0}(7 f+1)(2 f+1)^{5 / 2}.
   ```

2. `BirchMurnaghan3rd`:

   ```math
   B(V) = B_{0}(2 f+1)^{5 / 2} \left\{ 1 + (3B_{0}^{\prime} - 5) f + \frac{ 27 }{ 2 }(B_{0}^{\prime} - 4) f^2 \right\}
   ```

3. `BirchMurnaghan4th`:

   ```math
   B(V) = \frac{1}{6} B_{0}(2 f+1)^{5 / 2}\left\{\left(99 H-693 B_{0}^{\prime}+1573\right) f^{3}\right.\left.+\left(27 H-108 B_{0}^{\prime}+105\right) f^{2}+6\left(3 B_{0}^{\prime}-5\right) f+6\right\}.
   ```

4. `PoirierTarantola2nd`:

   ```math
   B(V) = \frac{B_{0}}{x}(1-\ln x).
   ```

5. `PoirierTarantola3rd`:

   ```math
   B(V) = -\frac{B_{0}}{2 x}\left[\left(B_{0}^{\prime}+2\right) \ln x(\ln x-1)-2\right].
   ```

6. `PoirierTarantola4th`:

   ```math
   B(V) = -\frac{B_{0}}{6 x}\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{3} x-3\left(H+2 B_{0}^{\prime}+1\right) \ln ^{2} x\right.\left.-6\left(B_{0}^{\prime}+1\right) \ln x-6\right\}.
   ```

7. `Vinet`:

   ```math
   B(V) = -\frac{B_{0}}{2 \eta^{2}}\left[3 \eta(\eta-1)\left(B_{0}^{\prime}-1\right)+2(\eta-2)\right]\times \exp \left\{-\frac{3}{2}\left(B_{0}^{\prime}-1\right)(\eta-1)\right\}.
   ```

8. `AntonSchmidt`:

   ```math
   B(V) = \beta\left(\frac{V}{V_{0}}\right)^{n}\left[1+n \ln \frac{V}{V_{0}}\right].
   ```

## Public interfaces

```@docs
Murnaghan1st
BirchMurnaghan
BirchMurnaghan2nd
BirchMurnaghan3rd
BirchMurnaghan4th
PoirierTarantola
PoirierTarantola2nd
PoirierTarantola3rd
Vinet
EnergyEquation
PressureEquation
BulkModulusEquation
getparam
orderof
real
isreal
float
ustrip
```
