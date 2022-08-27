```@meta
CurrentModule = EquationsOfStateOfSolids.FiniteStrains
```

# Finite strains

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

```@docs
VolumeTo
VolumeFrom
```
