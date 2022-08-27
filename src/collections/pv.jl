using Unitful: ħ, me

const FERMI_GAS_CONSTANT = (3π^2)^(2 / 3) * ħ^2 / 5 / me

function (eos::PressureEquation{<:Murnaghan1st})(v)
    @unpack v₀, b₀, b′₀ = eos
    return b₀ / b′₀ * ((v₀ / v)^b′₀ - 1)
end
function (eos::PressureEquation{<:Murnaghan2nd})(v)
    @unpack v₀, b₀, b′₀, b″₀ = eos
    h = b′₀^2 - 2b₀ * b″₀
    r = (v₀ / v)^h
    return 2b₀ / b′₀ / (h / b′₀ * ((r + 1) / (r - 1)) - 1)
end
function (eos::PressureEquation{<:BirchMurnaghan2nd})(v)
    @unpack v₀, b₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    return 3b₀ * f * (2f + 1)^_2½
end
function (eos::PressureEquation{<:BirchMurnaghan3rd})(v)
    @unpack v₀, b₀, b′₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    return 3f / 2 * b₀ * (2f + 1)^_2½ * (2 + 3f * (b′₀ - 4))
end
function (eos::PressureEquation{<:BirchMurnaghan4th})(v)
    @unpack v₀, b₀, b′₀, b″₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return b₀ / 2 * (2f + 1)^_2½ * ((9h - 63b′₀ + 143) * f^2 + 9f * (b′₀ - 4) + 6)
end
function (eos::PressureEquation{<:PoirierTarantola2nd})(v)
    @unpack v₀, b₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    return -3b₀ * f * exp(-3f)
end
function (eos::PressureEquation{<:PoirierTarantola3rd})(v)
    @unpack v₀, b₀, b′₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    return -3b₀ / 2 * f * exp(-3f) * (3f * (b′₀ - 2) + 1)
end
function (eos::PressureEquation{<:PoirierTarantola4th})(v)
    @unpack v₀, b₀, b′₀, b″₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return -3b₀ / 2 * f * exp(-3f) * (3f^2 * (h + 3b′₀ + 3) + 3f * (b′₀ - 2) + 2)
end
function (eos::PressureEquation{<:Vinet})(v)
    @unpack v₀, b₀, b′₀ = eos
    x, y = (v / v₀)^_⅓, 3 / 2 * (b′₀ - 1)
    return 3b₀ / x^2 * (1 - x) * exp(y * (1 - x))
end
function (eos::PressureEquation{<:AntonSchmidt})(v)
    @unpack v₀, b₀, b′₀ = eos
    x, n = v / v₀, -b′₀ / 2
    return -b₀ * x^n * log(x)
end
function (eos::PressureEquation{<:Holzapfel{Z}})(v) where {Z}
    @unpack v₀, b₀, b′₀ = eos
    η = (v / v₀)^_⅓
    p0 = FERMI_GAS_CONSTANT * (Z / v₀)^(5 / 3)
    c0 = -log(NoUnits(3b₀ / p0))
    c2 = 3 / 2 * (b′₀ - 3) - c0
    return 3b₀ * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
end
