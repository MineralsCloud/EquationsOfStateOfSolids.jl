module Volume

using UnPack: @unpack

using ..Collections

export volumeof

function volumeof(eos::PressureEoss{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
end
