using PolynomialRoots: roots
using Polynomials: fit as polyfit, derivative, coeffs

using ..EquationsOfStateOfSolids: FiniteStrainParameters, orderof
using ..FiniteStrains: FiniteStrain, VolumeTo, VolumeFrom, Dⁿᵥf, straintype

export linfit

struct LinearFitting <: FittingMethod end

struct CriterionNotMet
    msg::String
end

function linfit(volumes, energies, initial_params; kwargs...)
    return fit(volumes, energies, initial_params, LinearFitting(); kwargs...)
end

"""
    fit(volumes, energies, initial_params::FiniteStrainParameters, LinearFitting(); kwargs...)

Fit an equation of state ``E(V)`` using linear algorithms.

# Arguments
- `maxiter::Integer=1000`: .
- `conv_thr::AbstractFloat=1e-12`: .
- `root_thr::AbstractFloat=1e-20`: .
- `verbose::Bool=false`: .

!!! note
    If you want to fit with `BigFloat` data, you need to install
    [`GenericSVD.jl`](https://github.com/JuliaLinearAlgebra/GenericSVD.jl) and `using GenericSVD`
    before fittting!
"""
function fit(
    volumes,
    energies,
    initial_params::FiniteStrainParameters,
    ::LinearFitting;
    maxiter=1000,
    conv_thr=1e-12,
    root_thr=1e-20,
    verbose=false,
)
    deg, S = orderof(initial_params), straintype(initial_params)
    v0 = iszero(initial_params.v0) ? volumes[findmin(energies)[2]] : initial_params.v0  # Initial v0
    uv, ue = unit(v0), unit(energies[1])
    v0 = ustrip(v0)
    volumes = collect(map(x -> ustrip(uv, x), volumes))  # `parent` is needed to unwrap `DimArray`
    energies = collect(map(x -> ustrip(ue, x), energies))
    for i in 1:maxiter  # Self consistent loop
        strains = map(VolumeTo{S}(v0), volumes)
        if !(isreal(strains) && isreal(energies))
            throw(DomainError("the strains or the energies are complex!"))
        end
        poly = polyfit(real(strains), real(energies), deg)
        f0, e0 = min_of_min(poly, root_thr)
        v0_prev, v0 = v0, VolumeFrom{S}(v0)(f0)  # Record v0 to v0_prev, then update v0
        if abs((v0_prev - v0) / v0_prev) <= conv_thr
            if verbose
                @info "convergence reached after $i steps!"
            end
            fᵥ = map(deg -> Dⁿᵥf(S(), deg, v0)(v0), 1:4)
            e_f = map(deg -> derivative(poly, deg)(f0), 1:4)
            b0, b′0, b″0 = Dₚb(fᵥ, e_f)
            param = update(
                initial_params;
                v0=v0 * uv,
                b0=b0(v0) * ue / uv,
                b′0=b′0(v0),
                b″0=b″0(v0) * uv / ue,
                e0=e0 * ue,
            )
            return param
        end
    end
    throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
end

function islocalmin(x, y)  # `x` & `y` are both real
    y″ₓ = real(derivative(y, 2)(x))  # Must be real
    return y″ₓ > zero(y″ₓ)  # If 2nd derivative at x > 0, (x, y(x)) is a local minimum
end

function localmin(y, root_thr=1e-20)  # `y` is a polynomial (could be complex)
    y′ = derivative(y, 1)
    pool = roots(coeffs(y′); polish=true, epsilon=root_thr)
    real_roots = real(filter(isreal, pool))  # Complex volumes are meaningless
    if isempty(real_roots)
        throw(CriterionNotMet("no real extrema found! Consider changing `root_thr`!"))  # For some polynomials, could be all complex
    else
        localminima = filter(x -> islocalmin(x, y), real_roots)
        if isempty(localminima)
            throw(CriterionNotMet("no local minima found!"))
        else
            return localminima
        end
    end
end

# https://stackoverflow.com/a/21367608/3260253
function min_of_min(y, root_thr=1e-20)  # Find the minimum of the local minima
    localminima = localmin(y, root_thr)
    y0, i = findmin(map(y, localminima))  # `y0` must be real, or `findmap` will error
    x0 = localminima[i]
    return x0, y0
end

# See Eq. (55) - (57) in Ref. 1.
function Dₚb(fᵥ, e_f)  # Bulk modulus & its derivatives
    e″ᵥ = D²ᵥe(fᵥ, e_f)
    e‴ᵥ = D³ᵥe(fᵥ, e_f)
    b0 = v -> v * e″ᵥ
    b′0 = v -> -v * e‴ᵥ / e″ᵥ - 1
    b″0 = v -> (v * (D⁴ᵥe(fᵥ, e_f) * e″ᵥ - e‴ᵥ^2) + e‴ᵥ * e″ᵥ) / e″ᵥ^3
    return b0, b′0, b″0  # 3 lazy functions
end

function update(x::FiniteStrainParameters; kwargs...)
    patch = (; (f => kwargs[f] for f in propertynames(x))...)
    return setproperties(x, patch)
end

# Energy-volume derivatives, see Eq. (50) - (53) in Ref. 1.
D¹ᵥe(fᵥ, e_f) = e_f[1] * fᵥ[1]
D²ᵥe(fᵥ, e_f) = e_f[2] * fᵥ[1]^2 + e_f[1] * fᵥ[2]
D³ᵥe(fᵥ, e_f) = e_f[3] * fᵥ[1]^3 + 3fᵥ[1] * fᵥ[2] * e_f[2] + e_f[1] * fᵥ[3]
function D⁴ᵥe(fᵥ, e_f)
    return e_f[4] * fᵥ[1]^4 +
           6fᵥ[1]^2 * fᵥ[2] * e_f[3] +
           (4fᵥ[1] * fᵥ[3] + 3fᵥ[3]^2) * e_f[2] +
           e_f[1] * fᵥ[4]
end
