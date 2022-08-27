function (eos::EnergyEquation{<:Murnaghan1st})(v)
    @unpack v₀, b₀, b′₀, e₀ = eos
    x, y = b′₀ - 1, (v₀ / v)^b′₀
    return e₀ + b₀ / b′₀ * v * (y / x + 1) - v₀ * b₀ / x
end
function (eos::EnergyEquation{<:BirchMurnaghan2nd})(v)
    @unpack v₀, b₀, e₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    return e₀ + 9b₀ * v₀ * f^2 / 2
end
function (eos::EnergyEquation{<:BirchMurnaghan3rd})(v)
    @unpack v₀, b₀, b′₀, e₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    return e₀ + 9b₀ * v₀ * f^2 / 2 * (1 + f * (b′₀ - 4))
end
function (eos::EnergyEquation{<:BirchMurnaghan4th})(v)
    @unpack v₀, b₀, b′₀, b″₀, e₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return e₀ + 3b₀ * v₀ / 8 * f^2 * ((9h - 63b′₀ + 143) * f^2 + 12f * (b′₀ - 4) + 12)
end
function (eos::EnergyEquation{<:PoirierTarantola2nd})(v)
    @unpack v₀, b₀, e₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    return e₀ + 9b₀ * v₀ * f^2 / 2
end
function (eos::EnergyEquation{<:PoirierTarantola3rd})(v)
    @unpack v₀, b₀, b′₀, e₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    return e₀ + 9b₀ * v₀ * f^2 / 2 * ((2 - b′₀) * f + 1)
end
function (eos::EnergyEquation{<:PoirierTarantola4th})(v)
    @unpack v₀, b₀, b′₀, b″₀, e₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return e₀ + 9b₀ * v₀ * f^2 * (3f^2 * (h + 3b′₀ + 3) + 4f * (b′₀ + 2) + 4)
end
function (eos::EnergyEquation{<:Vinet})(v)
    @unpack v₀, b₀, b′₀, e₀ = eos
    x, y = 1 - (v / v₀)^_⅓, 3 / 2 * (b′₀ - 1)
    return e₀ + 9b₀ * v₀ / y^2 * (1 + (x * y - 1) * exp(x * y))
end
function (eos::EnergyEquation{<:AntonSchmidt})(v)
    @unpack v₀, b₀, b′₀, e∞ = eos
    n₊₁ = -b′₀ / 2 + 1
    x = v / v₀
    return e∞ + b₀ * v₀ / n₊₁ * x^n₊₁ * (log(x) - 1 / n₊₁)
end
