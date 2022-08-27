(eos::BulkModulusEquation{<:Murnaghan1st})(v) = eos.b₀ + PressureEquation(eos)(v)
function (eos::BulkModulusEquation{<:BirchMurnaghan2nd})(v)
    @unpack v₀, b₀ = eos
    f = EulerianStrainFromVolume(v₀)(v)
    return b₀ * (7f + 1) * (2f + 1)^_2½
end
# The formula in the Gibbs2 paper is wrong! See the EosFit7c paper!
function (eos::BulkModulusEquation{<:BirchMurnaghan3rd})(v)
    v₀, b₀, b′₀ = unpack(eos)
    f = EulerianStrainFromVolume(v₀)(v)
    return b₀ * (2f + 1)^_2½ * (27f^2 / 2 * (b′₀ - 4) + (3b′₀ - 5) * f + 1)
end
function (eos::BulkModulusEquation{<:BirchMurnaghan4th})(v)
    v₀, b₀, b′₀, b″₀ = unpack(eos)
    f = EulerianStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return b₀ / 6 *
           (2f + 1)^_2½ *
           ((99h - 693b′₀ + 1573) * f^3 + (27h - 108b′₀ + 105) * f^2 + 6f * (3b′₀ - 5) + 6)
end
function (eos::BulkModulusEquation{<:PoirierTarantola2nd})(v)
    @unpack v₀, b₀ = eos
    f = NaturalStrainFromVolume(v₀)(v)
    return b₀ * (1 - 3f) * exp(-3f)
end
function (eos::BulkModulusEquation{<:PoirierTarantola3rd})(v)
    v₀, b₀, b′₀ = unpack(eos)
    f = NaturalStrainFromVolume(v₀)(v)
    return -b₀ / 2 * exp(-3f) * (9f^2 * (b′₀ - 2) - 6f * (b′₀ + 1) - 2)
end
function (eos::BulkModulusEquation{<:PoirierTarantola4th})(v)
    v₀, b₀, b′₀, b″₀ = unpack(eos)
    f = NaturalStrainFromVolume(v₀)(v)
    h = b″₀ * b₀ + b′₀^2
    return -b₀ / 2 *
           exp(-3f) *
           (9f^3 * (h + 3b′₀ + 3) - 9f^2 * (h + 2b′₀ + 1) - 6f * (b′₀ + 1) - 2)
end
function (eos::BulkModulusEquation{<:Vinet})(v)
    v₀, b₀, b′₀ = unpack(eos)
    x, ξ = (v / v₀)^_⅓, 3 / 2 * (b′₀ - 1)
    return -b₀ / (2 * x^2) * (3x * (x - 1) * (b′₀ - 1) + 2 * (x - 2)) * exp(-ξ * (x - 1))
end
function (eos::BulkModulusEquation{<:AntonSchmidt})(v)
    v₀, b₀, b′₀ = unpack(eos)
    x, n = v / v₀, -b′₀ / 2
    return b₀ * x^n * (1 + n * log(x))
end
