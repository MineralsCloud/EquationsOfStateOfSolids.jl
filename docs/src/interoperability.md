# How to use `EquationsOfStateOfSolids` in Python?

It may be attempting for [Pythonistas](https://en.wiktionary.org/wiki/Pythonista)
to use this package in Python, without
writing too much code. Luckily, Julia provides such a feature.

1. First, install [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl), following their [instructions](https://github.com/JuliaPy/PyCall.jl/blob/master/README.md). Notice on macOS, that if you want to install Python from [`pyenv`](https://github.com/pyenv/pyenv), you may need to run

   ```shell
   env PYTHON_CONFIGURE_OPTS="--enable-framework" pyenv install 3.6.9
   ```

   in terminal to install your `python`, or else Julia will throw an

   ```julia
   ImportError: No module named site
   ```

   See [this issue](https://github.com/JuliaPy/PyCall.jl/issues/122) and [another issue](https://github.com/JuliaPy/PyCall.jl/issues/597) for details.

2. Install [`PyJulia`](https://pyjulia.readthedocs.io/en/stable/index.html) in Python. Please see [its official tutorial](https://pyjulia.readthedocs.io/en/stable/installation.html#step-2-install-pyjulia) for instructions.

3. Open a (an) Python (IPython) session, start playing!

   ```python
   In [1]: from julia import Unitful

   In [2]: from julia.EquationsOfStateOfSolids.Collections import *

   In [3]: from julia.EquationsOfStateOfSolids.Fitting import *

   In [4]: Murnaghan(1, 2, 3.0, 4)
   Out[4]: <PyCall.jlwrap EquationsOfStateOfSolids.Collections.Murnaghan{Float64}(1.0, 2.0, 3.0, 4.0)>

   In [5]: result = nonlinfit(
      ...:     PressureEquation(BirchMurnaghan3rd(1, 2, 3.0, 0)),
      ...:     [1, 2, 3, 4, 5],
      ...:     [5, 6, 9, 8, 7],
      ...: )

   In [6]: result.v0, result.b0, result.bp0
   Out[6]: (1.1024687826913997, 29.308616965851673, 12.689089874230556)

   In [7]: from julia import Main

   In [8]: volumes = Main.eval("data[:, 1] .* UnitfulAtomic.bohr^3")

   In [9]: energies = Main.eval("data[:, 2] .* UnitfulAtomic.Ry")
   ```

   where `data` is copied from Julia:

   ```python
   In [1]: data = Main.eval("""
      ...:    [
      ...:        159.9086 -323.4078898
                      ⋮          ⋮
      ...:        319.8173 -323.4105393
      ...:    ]
      ...:    """
      ...: )
   ```
