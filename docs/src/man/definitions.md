# Definitions and conventions

```@contents
Pages = ["definitions.md"]
Depth = 2
```

## Finite strains

This module contains some methods to calculate several finite strains.
The following formulae are from
[the `Gibbs2` paper](https://www.sciencedirect.com/science/article/pii/S0010465511001470) Table 3.

- Eulerian strain:

  ```math
  f = \frac{ 1 }{ 2 } \bigg( \Big \frac{ V_0 }{ V } \Big^{\frac{ 2 }{ 3 }} - 1 \bigg)
  ```

- Lagrangian strain:

  ```math
  f = \frac{ 1 }{ 2 } \bigg( \Big( \frac{ V }{ V_0 } \Big^{\frac{ 2 }{ 3 }} - 1 \bigg)
  ```

- Natural (Hencky) strain:

  ```math
  f = \frac{ 1 }{ 3 } \ln \Big \frac{ V }{ V_0 } \Big
  ```

- Infinitesimal strain:

  ```math
  f = 1 - \Big \frac{ V_0 }{ V } \Big^{\frac{ 1 }{ 3 }}
  ```

## Energy equations of state

The ``E(V)`` relation of equations of state are listed as below:

1. `Murnaghan`:

2. `Birch--Murnaghan 2nd`:

   ```math
   ```

3. `Birch--Murnaghan 3rd`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

   where ``x = V / V_0``, and
   ``f = \frac{ 1 }{ 2 } \bigg[ \bigg( \frac{ V_0 }{ V } \bigg)^{2/3} - 1 \bigg]``.

4. `Birch--Murnaghan 4th`:

   ```math
   E(V) = E_{0}+\frac{3}{8} V_{0} B_{0} f^{2}\left[\left(9 H-63 B_{0}^{\prime}+143\right) f^{2}+12\left(B_{0}^{\prime}-4\right) f+12\right].
   ```

   where ``H = B_0 B_0'' + (B_0')^2``.

5. `Poirier--Tarantola 2nd`:

   ```math
   E(V) = E_{0}+\frac{1}{2} B_{0} V_{0} \ln ^{2} x.
   ```

6. `Poirier--Tarantola 3rd`:

   ```math
   E(V) = E_{0}+\frac{1}{6} B_{0} V_{0} \ln ^{2} x\left[\left(B_{0}^{\prime}+2\right) \ln x+3\right].
   ```

7. `Poirier--Tarantola 4th`:

   ```math
   E(V) = E_{0}+\frac{1}{24} B_{0} V_{0} \ln ^{2} x\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{2} x\right. \left.+4\left(B_{0}^{\prime}+2\right) \ln x+12\right\}.
   ```

   where ``H = B_0 B_0'' + (B_0')^2``.

8. `Vinet`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

9. `Anton--Schmidt`:

   ```math
   E(V)=\frac{\beta V_{0}}{n+1}\left(\frac{V}{V_{0}}\right)^{n+1}\left[\ln \left(\frac{V}{V_{0}}\right)-\frac{1}{n+1}\right]+E_{\infty}.
   ```

## Pressure equations of state

The ``P(V)`` relation of equations of state are listed as below:

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

## Bulk modulus equations of state

The ``B(V)`` relation of equations of state are listed as below:

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

## Linear fitting

The linear fitting

## Nonlinear fitting

> The equations of state depend nonlinearly on a collection of parameters,
> $E_0$, $V_0$, $B_0$, $B_0'$, ..., that represent physical properties of the
> solid at equilibrium and can, in principle, be obtained experimentally by
> independent methods. The use of a given analytical EOS may have significant
> influence on the results obtained, particularly because the parameters are far
> from being independent. The number of parameters has to be considered in
> comparing the goodness of fit of functional forms with different analytical
> flexibility. The possibility of using too many parameters, beyond what is
> physically justified by the information contained in the experimental data, is
> a serious aspect that deserves consideration.[^1]

In [`EquationsOfStateOfSolids`](https://github.com/MineralsCloud/EquationsOfStateOfSolids.jl),
the nonlinear fitting is currently implemented by
[`LsqFit`](https://github.com/JuliaNLSolvers/LsqFit.jl), a small library that
provides basic least-squares fitting in pure Julia. It only utilizes the
_Levenberg–Marquardt algorithm_ for non-linear fitting. See its
[documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md)
for more information.
