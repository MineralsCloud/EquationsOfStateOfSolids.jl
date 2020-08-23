module Volume

using InteractiveUtils: subtypes
using PolynomialRoots: roots
using Roots: find_zero, AbstractBracketing, AbstractSecant
using UnPack: @unpack

using ..Collections:
    PressureEOS,
    EnergyEOS,
    EquationOfStateOfSolids,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    PoirierTarantola3rd,
    EulerianStrain,
    strain2volume

export findvolume

function _allsubtypes(t::Type, types = Type[])
    for s in subtypes(t)
        types = _allsubtypes(s, push!(types, s))
    end
    return types
end
_nonabstract(t::Type) = filter(!isabstracttype, _allsubtypes(t))  # `Roots.FalsePosition` is not concrete!

const ROOT_FINDING_ALGORITHMS = _nonabstract(AbstractUnivariateZeroMethod)

function findvolume(eos::PressureEOS{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function findvolume(eos::EnergyEOS{<:BirchMurnaghan2nd}, e)
    @unpack v0, b0, b′0, e0 = eos.param
    f = sqrt(2 / 9 * (e - e0) / b0 / v0)
    vs = map(strain2volume(EulerianStrain(), v0), [f, -f])
    return map(real, filter(isreal, vs))
end
function findvolume(eos::EnergyEOS{<:BirchMurnaghan3rd}, e; root_thr = 1e-20)
    @unpack v0, b0, b′0, e0 = eos.param
    # Constrcut ax^3 + bx^2 + d = 0
    b, d = 9 / 2 * b0 * v0, e0 - e
    a = b * (b′0 - 4)
    # Solve ax^3 + bx^2 + d = 0
    fs = roots([d, 0, b, a]; polish = true, epsilon = root_thr)
    vs = map(strain2volume(EulerianStrain(), v0), fs)
    return map(real, filter(isreal, vs))
end
function findvolume(
    eos::EquationOfStateOfSolids,
    y,
    method::Union{AbstractBracketing,AbstractSecant};
    volume_scale = (0.5, 1.5),
    maxiter = 40,
    verbose = false,
)
    v = find_zero(
        x -> eos(x) - y,
        _volume_scale(volume_scale, method) .* eos.param.v0,
        method;
        maxevals = maxiter,
        verbose = verbose,
    )
    if v < zero(v)
        @warn "the volume found is negative!"
    end
    return v
end
_volume_scale(volume_scale, ::AbstractBracketing) = extrema(volume_scale)
_volume_scale(volume_scale, ::AbstractSecant) = sum(extrema(volume_scale)) / 2

end
