(eos::BulkModulusEquation{<:Murnaghan1st})(v) =
    getparam(eos).b0 + PressureEquation(getparam(eos))(v)
function (eos::BulkModulusEquation{<:BirchMurnaghan2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    return b0 * (7f + 1) * (2f + 1)^_2½
end
function (eos::BulkModulusEquation{<:BirchMurnaghan3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    return b0 / 2 * (2f + 1)^_2½ * ((27f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (eos::BulkModulusEquation{<:BirchMurnaghan4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = VolumeToEulerianStrain(v0)(v)
    h = b″0 * b0 + b′0^2
    return b0 / 6 *
           (2f + 1)^_2½ *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function (eos::BulkModulusEquation{<:PoirierTarantola2nd})(v)
    @unpack v0, b0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    return b0 * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEquation{<:PoirierTarantola3rd})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    return -b0 / 2 * exp(-3f) * (9f^2 * (b′0 - 2) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEquation{<:PoirierTarantola4th})(v)
    @unpack v0, b0, b′0, b″0 = getparam(eos)
    f = VolumeToNaturalStrain(v0)(v)
    h = b″0 * b0 + b′0^2
    return -b0 / 2 *
           exp(-3f) *
           (9f^3 * (h + 3b′0 + 3) - 9f^2 * (h + 2b′0 + 1) - 6f * (b′0 + 1) - 2)
end
function (eos::BulkModulusEquation{<:Vinet})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, ξ = (v / v0)^_⅓, 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end
function (f::BulkModulusEquation{<:AntonSchmidt})(v)
    @unpack v0, b0, b′0 = getparam(eos)
    x, n = v / v0, -b′0 / 2
    return b0 * x^n * (1 + n * log(x))
end
