module Volume

using Compat: filter
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

using EquationsOfStateOfSolids: _ispositive
using ..Collections:
    PressureEos,
    EnergyEos,
    EquationOfStateOfSolids,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    PoirierTarantola3rd,
    EulerianStrain,
    strain2volume,
    getparam

export findvolume, mustfindvolume

function findvolume(eos::PressureEos{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function findvolume(eos::EnergyEos{<:BirchMurnaghan2nd}, e)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    f = sqrt(2 / 9 * (e - e0) / b0 / v0)
    vs = map(strain2volume(EulerianStrain(), v0), (f, -f))
    return map(real, filter(isreal, vs))
end
function findvolume(eos::EnergyEos{<:BirchMurnaghan3rd}, e)
    @unpack v0, b0, b′0, e0 = getparam(eos)
    # Constrcut ax^3 + bx^2 + d = 0, see https://zh.wikipedia.org/wiki/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
    a = b′0 - 4
    r = 1 / 3a  # b = 1
    d = (e0 - e) / (9 / 2 * b0 * v0)
    p, q = -r^3 - d / 2a, -r^2
    Δ = p^2 + q^3
    fs = -r .+ if Δ < 0
        SIN, COS = sincos(acos(p / abs(r)^3) / 3)
        (2COS, -COS - √3 * SIN, -COS + √3 * SIN) .* abs(r)  # Verified
    elseif Δ == 0
        if p == q == 0
            (0,)  # 3 reals are equal
        else  # p == -q != 0
            2cbrt(p), -cbrt(p)  # 2 real roots are equal, leaving 2 solutions
        end
    else  # Δ > 0
        (cbrt(p + √Δ) + cbrt(p - √Δ),)  # Only 1 real solution
    end  # solutions are strains
    vs = map(strain2volume(EulerianStrain(), v0), fs)
    return filter(_ispositive, map(real, filter(isreal, vs)))
end
function findvolume(
    eos::EquationOfStateOfSolids,
    y,
    method::Union{AbstractBracketing,AbstractSecant};
    vscale = (0.5, 1.5),
    maxiter = 40,
    verbose = false,
)
    v0 = _vscale(vscale, method) .* getparam(eos).v0  # v0 can be negative
    @assert _ispositive(minimum(v0))  # No negative volume
    v = find_zero(x -> eos(x) - y, v0, method; maxevals = maxiter, verbose = verbose)
    if !_ispositive(v)
        @warn "the volume found is negative!"
    end
    return v
end
_vscale(vscale, ::AbstractBracketing) = extrema(vscale)
_vscale(vscale, ::AbstractSecant) = sum(extrema(vscale)) / 2

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
    error("no volume found!")
end

end
