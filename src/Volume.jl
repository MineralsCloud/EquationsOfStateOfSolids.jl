module Volume

using PolynomialRoots: roots
using UnPack: @unpack

using ..Collections

export volumeof, volumeof2, volumeof3

function volumeof(eos::PressureEoss{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function volumeof(eos::EnergyEoss{<:BirchMurnaghan3rd}, e)
    @unpack v0, b0, b′0, e0 = eos.param  # https://math.vanderbilt.edu/schectex/courses/cubic/
    # Constrcut ax^3 + bx^2 + d = 0
    b, d = 9 / 2 * b0 * v0, e0 - e
    a = b * (b′0 - 4)
    # Solve ax^3 + bx^2 + d = 0
    p = -b / 3 / a
    q = p^3 - d / 2 / a
    r = sqrt(Complex(q^2 - p^6))
    f = (q + r)^(1 / 3) + (q - r)^(1 / 3) + p  # The real solution
    fs = f .* (1, (-1 + √3im) / 2, (-1 - √3im) / 2)
    vs = filter(isreal, map(volume_from_strain(Eulerian(), v0), fs))
    return map(real, vs)
end
function volumeof3(eos::EnergyEoss{<:BirchMurnaghan3rd}, e)
    @unpack v0, b0, b′0, e0 = eos.param
    # Constrcut ax^3 + bx^2 + d = 0
    b, d = 9 / 2 * b0 * v0, e0 - e
    a = b * (b′0 - 4)
    # Solve ax^3 + bx^2 + d = 0
    fs = roots([d, 0, b, a], polish = true, epsilon = 1e-20)
    vs = map(volume_from_strain(Eulerian(), v0), fs)
    return map(real, filter(isreal, vs))
end

end
