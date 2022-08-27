using Unitful: ħ, me

const FERMI_GAS_CONSTANT = (3π^2)^(2 / 3) * ħ^2 / 5 / me

function (eos::PressureEquation{<:Murnaghan1st})(v)
    @unpack v0, b0, b′0 = eos
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function (eos::PressureEquation{<:Murnaghan2nd})(v)
    @unpack v0, b0, b′0, b″0 = eos
    h = b′0^2 - 2b0 * b″0
    r = (v0 / v)^h
    return 2b0 / b′0 / (h / b′0 * ((r + 1) / (r - 1)) - 1)
end
function (eos::PressureEquation{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = eos
    f = EulerianStrainFromVolume(v0)(v)
    return 3b0 * f * (2f + 1)^_2½
end
function (eos::PressureEquation{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = eos
    f = EulerianStrainFromVolume(v0)(v)
    return 3f / 2 * b0 * (2f + 1)^_2½ * (2 + 3f * (b′0 - 4))
end
function (eos::PressureEquation{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = eos
    f = EulerianStrainFromVolume(v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 2 * (2f + 1)^_2½ * ((9h - 63b′0 + 143) * f^2 + 9f * (b′0 - 4) + 6)
end
function (eos::PressureEquation{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = eos
    f = NaturalStrainFromVolume(v0)(v)
    return -3b0 * f * exp(-3f)
end
function (eos::PressureEquation{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = eos
    f = NaturalStrainFromVolume(v0)(v)
    return -3b0 / 2 * f * exp(-3f) * (3f * (b′0 - 2) + 1)
end
function (eos::PressureEquation{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = eos
    f = NaturalStrainFromVolume(v0)(v)
    h = b″0 * b0 + b′0^2
    return -3b0 / 2 * f * exp(-3f) * (3f^2 * (h + 3b′0 + 3) + 3f * (b′0 - 2) + 2)
end
function (eos::PressureEquation{<:Vinet})(v)
    @unpack v0, b0, b′0 = eos
    x, y = (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end
function (eos::PressureEquation{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = eos
    x, n = v / v0, -b′0 / 2
    return -b0 * x^n * log(x)
end
function (eos::PressureEquation{<:Holzapfel{Z}})(v) where {Z}
    @unpack v0, b0, b′0 = eos
    η = (v / v0)^_⅓
    p0 = FERMI_GAS_CONSTANT * (Z / v0)^(5 / 3)
    c0 = -log(NoUnits(3b0 / p0))
    c2 = 3 / 2 * (b′0 - 3) - c0
    return 3b0 * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
end
