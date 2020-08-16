module Volume

using PolynomialRoots: roots
using Roots:
    find_zero,
    Halley,
    Schroder,
    Bisection,
    BisectionExact,
    FalsePosition,
    A42,
    AlefeldPotraShi,
    Esser,
    King,
    KumarSinghAkanksha,
    Order0,
    Order16,
    Order1B,
    Order2,
    Order2B,
    Order5,
    Order8,
    Secant,
    Steffensen,
    Thukral16,
    Thukral8,
    Brent,
    Newton
using UnPack: @unpack

using ..Collections

export findvolume

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
function findvolume(eos::EquationOfStateOfSolids, y, v0, method)
    v = _findvolume(eos, y, v0, method)
    if v < zero(v)
        @error "the volume found is negative!"
    end
    return v
end # function findvolume
function findvolume(eos::EquationOfStateOfSolids, y, v0; silent = false)
    for T in [
        Bisection,
        BisectionExact,
        FalsePosition,
        A42,
        AlefeldPotraShi,
        Brent,
        Halley,
        Schroder,
        Newton,
        Esser,
        King,
        KumarSinghAkanksha,
        Order0,
        Order16,
        Order1B,
        Order2,
        Order2B,
        Order5,
        Order8,
        Secant,
        Steffensen,
        Thukral16,
        Thukral8,
    ]
        silent || @info "using method `$T`..."
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return findvolume(eos, y, v0, T())
        catch e
            silent || @info "method `$T` failed because of $e."
            continue
        end
    end
end # function findvolume
_findvolume(
    eos,
    y,
    v0,
    method::Union{Bisection,BisectionExact,FalsePosition,A42,AlefeldPotraShi},
) = find_zero(v -> eos(v) - y, (minimum(v0), maximum(v0)), method)
_findvolume(
    eos,
    y,
    v0,
    method::Union{
        Brent,
        Halley,
        Schroder,
        Newton,
        Esser,
        King,
        KumarSinghAkanksha,
        Order0,
        Order16,
        Order1B,
        Order2,
        Order2B,
        Order5,
        Order8,
        Secant,
        Steffensen,
        Thukral16,
        Thukral8,
    },
) = find_zero(v -> eos(v) - y, (minimum(v0) + maximum(v0)) / 2, method)

end
