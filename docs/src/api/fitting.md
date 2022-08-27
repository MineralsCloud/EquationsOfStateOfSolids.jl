```@meta
CurrentModule = EquationsOfStateOfSolids.Fitting
```

# Fitting

```@contents
Pages = ["fitting.md"]
```

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

## Linear fitting

The linear fitting

## Usage

```@repl
using EquationsOfStateOfSolids
using EquationsOfStateOfSolids.Fitting
volumes = [
    25.987454833,
    26.9045702104,
    27.8430241908,
    28.8029649591,
    29.7848370694,
    30.7887887064,
    31.814968055,
    32.8638196693,
    33.9353435494,
    35.0299842495,
    36.1477417695,
    37.2892088485,
    38.4543854865,
    39.6437162376,
    40.857201102,
    42.095136449,
    43.3579668329,
    44.6456922537,
    45.9587572656,
    47.2973100535,
    48.6614988019,
    50.0517680652,
    51.4682660281,
    52.9112890601,
    54.3808371612,
    55.8775030703,
    57.4014349722,
    58.9526328669,
];
energies = [
    -7.63622156576,
    -8.16831294894,
    -8.63871612686,
    -9.05181213218,
    -9.41170988374,
    -9.72238224345,
    -9.98744832526,
    -10.210309552,
    -10.3943401353,
    -10.5427238068,
    -10.6584266073,
    -10.7442240979,
    -10.8027285713,
    -10.8363890521,
    -10.8474912964,
    -10.838157792,
    -10.8103477586,
    -10.7659387815,
    -10.7066179666,
    -10.6339907853,
    -10.5495538639,
    -10.4546677714,
    -10.3506386542,
    -10.2386366017,
    -10.1197772808,
    -9.99504030111,
    -9.86535084973,
    -9.73155247952,
];
nonlinfit(volumes, energies, BirchMurnaghan3rd(40, 0.5, 4))
nonlinfit(volumes, energies, Murnaghan(41, 0.5, 4))
nonlinfit(volumes, energies, PoirierTarantola3rd(41, 0.5, 4))
nonlinfit(volumes, energies, Vinet(41, 0.5, 4))
```

Then 4 different equations of state will be fitted.

They just work as well with units:

```@repl
using Unitful
volumes = volumes * u"angstrom^3"
energies = energies * u"eV"
nonlinfit(volumes, energies, BirchMurnaghan3rd(40u"angstrom^3", 1u"GPa", 4))
```

## Public interfaces

```@docs
fit
linfit
nonlinfit
```

## References

1. [A. Otero-De-La-Roza, V. Luaña, _Comput. Phys. Commun._ **182**, 1708–1720 (2011).](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
