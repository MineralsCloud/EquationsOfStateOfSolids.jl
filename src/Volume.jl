module Volume

using InteractiveUtils: subtypes
using PolynomialRoots: roots
using Roots:
    find_zero,
    AbstractBracketing,
    AbstractSecant,
    Bisection,
    BisectionExact,
    FalsePosition,
    A42,
    AlefeldPotraShi,
    Brent,
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
    Thukral8
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

export findvolume, mustfindvolume

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
    v0 = eos.param.v0  # v0 can be negative
    @assert _ispositive(volume_scale * v0)  # No negative volume
    v = find_zero(
        x -> eos(x) - y,
        _volume_scale(volume_scale, method) .* v0,
        method;
        maxevals = maxiter,
        verbose = verbose,
    )
    if !_ispositive(v)
        @warn "the volume found is negative!"
    end
    return v
end
_volume_scale(volume_scale, ::AbstractBracketing) = extrema(volume_scale)
_volume_scale(volume_scale, ::AbstractSecant) = sum(extrema(volume_scale)) / 2

function mustfindvolume(eos::EquationOfStateOfSolids, y; verbose = false, kwargs...)
    for T in [
        Bisection,
        BisectionExact,
        FalsePosition,
        A42,
        AlefeldPotraShi,
        Brent,
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
        if verbose
            @info "using method `$T`..."
        end
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return findvolume(eos, y, T(); verbose = verbose, kwargs...)
        catch e
            if verbose
                @info "method `$T` failed because of `$e`."
            end
            continue
        end
    end
end

_ispositive(x) = x > zero(x)

end
