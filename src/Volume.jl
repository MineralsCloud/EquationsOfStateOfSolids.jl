module Volume

using PolynomialRoots: roots
using Roots: find_zero, AbstractBracketing, AbstractUnivariateZeroMethod
using Statistics: mean
using UnPack: @unpack

using ..Collections:
    PressureEoss,
    EnergyEoss,
    EquationOfStateOfSolids,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd

export findvolume

function _allsubtypes(t::Type, types = Type[])
    for s in subtypes(t)
        types = _allsubtypes(s, push!(types, s))
    end
    return types
end
_nonabstract(t::Type) = filter(!isabstracttype, _allsubtypes(t))

const ROOT_FINDING_ALGORITHMS = _nonabstract(AbstractUnivariateZeroMethod)

function findvolume(eos::PressureEoss{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function findvolume(eos::EnergyEoss{<:BirchMurnaghan2nd}, e)
    @unpack v0, b0, b′0, e0 = eos.param
    f = sqrt(2 / 9 * (e - e0) / b0 / v0)
    vs = map(volume_from_strain(Eulerian(), v0), [f, -f])
    return map(real, filter(isreal, vs))
end
function findvolume(eos::EnergyEoss{<:BirchMurnaghan3rd}, e; epsilon = 1e-20)
    @unpack v0, b0, b′0, e0 = eos.param
    # Constrcut ax^3 + bx^2 + d = 0
    b, d = 9 / 2 * b0 * v0, e0 - e
    a = b * (b′0 - 4)
    # Solve ax^3 + bx^2 + d = 0
    fs = roots([d, 0, b, a]; polish = true, epsilon = epsilon)
    vs = map(volume_from_strain(Eulerian(), v0), fs)
    return map(real, filter(isreal, vs))
end
function findvolume(
    eos::EquationOfStateOfSolids,
    y,
    v0,
    method::AbstractUnivariateZeroMethod,
)
    v = _findvolume(eos, y, v0, method)
    if v < zero(v)
        @error "the volume found is negative!"
    end
    return v
end # function findvolume
function findvolume(eos::EquationOfStateOfSolids, y, v0; silent = false)
    for T in ROOT_FINDING_ALGORITHMS
        silent || @info "using method `$T`..."
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return findvolume(eos, y, v0, T())
        catch e
            silent || @info "method `$T` failed because of $e."
        end
    end
end # function findvolume
_findvolume(eos, y, v0, method::AbstractBracketing) =
    find_zero(v -> eos(v) - y, extrema(v0), method)
# The rest of root-finding methods
_findvolume(eos, y, v0, method::AbstractUnivariateZeroMethod) =
    find_zero(v -> eos(v) - y, mean(v0), method)

end
