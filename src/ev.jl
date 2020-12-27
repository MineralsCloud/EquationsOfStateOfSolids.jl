function (eos::EnergyEquation{<:Murnaghan1st})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    x, y = b′0 - 1, (v0 / v)^b′0
    return e0 + b0 / b′0 * v * (y / x + 1) - v0 * b0 / x
end
function (eos::EnergyEquation{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0, e0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
function (eos::EnergyEquation{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * (1 + f * (b′0 - 4))
end
function (eos::EnergyEquation{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 3b0 * v0 / 8 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12f * (b′0 - 4) + 12)
end
function (eos::EnergyEquation{<:PoirierTarantola2nd})(v)
    @unpack v0, b0, e0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2
end
"""
    (eos::EnergyEquation{<:PoirierTarantola3rd})(v)

Evaluate a Poirier--Tarantola energy EOS on volume `v`.

```math
E(V) = E_0 + \\frac {B_0 V_0}{2}\\left[ \\ln \\left( \\frac{V_0}{V}\\right) \\right]^{2} + \\frac{B_0 V_0}{6} \\left[ \\ln \\left(\\frac{V_0}V{} \\right) \\right]^3 \\left(B_0^\\prime - 2 \\right)
```
"""
function (eos::EnergyEquation{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    return e0 + 9b0 * v0 * f^2 / 2 * ((2 - b′0) * f + 1)
end
function (eos::EnergyEquation{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0, e0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    h = b″0 * b0 + b′0^2
    return e0 + 9b0 * v0 * f^2 * (3f^2 * (h + 3b′0 + 3) + 4f * (b′0 + 2) + 4)
end
function (eos::EnergyEquation{<:Vinet})(v)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    x, y = 1 - (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end
function (f::EnergyEquation{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0, e∞ = getparam(eos)
    n₊₁ = -b′0 / 2 + 1
    x = v / v0
    return e∞ + b0 * v0 / n₊₁ * x^n₊₁ * (log(x) - 1 / n₊₁)
end
