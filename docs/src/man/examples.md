# Usage

## Construct a `EquationOfStateOfSolidsParameters` instance

We will use `BirchMurnaghan3rd` as an example.

A `BirchMurnaghan3rd` can be constructed from scratch, as shown above. It can
also be constructed from an existing `BirchMurnaghan3rd`, with
[Setfield.jl](https://github.com/jw3126/Setfield.jl)
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
[Kaleido.jl](https://github.com/tkf/Kaleido.jl):

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
