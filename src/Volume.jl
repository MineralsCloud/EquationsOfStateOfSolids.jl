module Volume

using PolynomialRoots: roots
using UnPack: @unpack

using ..Collections

export volumeof

function volumeof(eos::PressureEoss{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function volumeof(eos::EnergyEoss{<:BirchMurnaghan3rd}, e, epsilon = 1e-20)
    @unpack v0, b0, b′0, e0 = eos.param
    # Constrcut ax^3 + bx^2 + d = 0
    b, d = 9 / 2 * b0 * v0, e0 - e
    a = b * (b′0 - 4)
    # Solve ax^3 + bx^2 + d = 0
    fs = roots([d, 0, b, a]; polish = true, epsilon = epsilon)
    vs = map(volume_from_strain(Eulerian(), v0), fs)
    return map(real, filter(isreal, vs))
end

end
