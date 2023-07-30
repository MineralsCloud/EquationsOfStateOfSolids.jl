using Polynomials: fit as polyfit

using ..EquationsOfStateOfSolids:
    FiniteStrainParameters, BirchMurnaghan, PoirierTarantola, orderof, straintype
using ..FiniteStrains:
    StrainFromVolume, VolumeFromStrain, VolumeFromStrain, DerivativeOfStrain

export linfit

struct LinearFitting <: FittingMethod end

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
    [GenericSVD.jl](https://github.com/JuliaLinearAlgebra/GenericSVD.jl) and `using GenericSVD`
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
    v₀ = iszero(initial_params.v0) ? volumes[argmin(energies)] : initial_params.v0  # Initial v0
    volumes = collect(volumes)
    energies = collect(energies)
    degree = orderof(initial_params)
    st = straintype(initial_params)
    for i in 1:maxiter  # Self consistent loop
        v₀_next, e₀ = swap(volumes, energies, degree, st, v₀, root_thr)
        v₀_prev, v₀ = v₀, v₀_next  # Record v0 to v0_prev, then update v0
        if abs((v₀_prev - v₀) / v₀_prev) <= conv_thr
            if verbose
                @info "convergence reached after $i steps!"
            end
            fᵥ = map(deg -> DerivativeOfStrain{typeof(st)}(deg)(v₀)(v), 1:4)
            e_f = map(deg -> derivative(poly, deg)(f0), 1:4)
            b0, b′0, b″0 = Dₚb(fᵥ, e_f)
            param = update(
                initial_params; v0=v₀, b0=b0(v₀), b′0=b′0(v₀), b″0=b″0(v₀), e0=e₀
            )
            return param
        end
    end
    throw(ConvergenceFailed("convergence not reached after $maxiter steps!"))
end

function swap(volumes, energies, degree, st, v₀, root_thr)
    strains = map(StrainFromVolume(st)(v₀), volumes)
    polynomial = polyfit(strains, energies, degree)
    f₀, e₀ = globalminimum(polynomial, root_thr)
    return VolumeFromStrain(st)(v₀)(f₀), e₀
end

# See Eq. (55) - (57) in Ref. 1.
function Dₚb(fᵥ, e_f)  # Bulk modulus & its derivatives
    e″ᵥ = DerivativeOfEnergy{2}()(fᵥ, e_f)
    e‴ᵥ = DerivativeOfEnergy{3}()(fᵥ, e_f)
    b0 = v -> v * e″ᵥ
    b′0 = v -> -v * e‴ᵥ / e″ᵥ - 1
    b″0 = v -> (v * (DerivativeOfEnergy{4}()(fᵥ, e_f) * e″ᵥ - e‴ᵥ^2) + e‴ᵥ * e″ᵥ) / e″ᵥ^3
    return b0, b′0, b″0  # 3 lazy functions
end

function update(x::FiniteStrainParameters; kwargs...)
    patch = (; (f => kwargs[f] for f in propertynames(x))...)
    return setproperties(x, patch)
end

# Energy-volume derivatives, see Eq. (50) - (53) in Ref. 1.
struct DerivativeOfEnergy{N} end
(::DerivativeOfEnergy{1})(fᵥ, e_f) = e_f[1] * fᵥ[1]
(::DerivativeOfEnergy{2})(fᵥ, e_f) = e_f[2] * fᵥ[1]^2 + e_f[1] * fᵥ[2]
function (::DerivativeOfEnergy{3})(fᵥ, e_f)
    return e_f[3] * fᵥ[1]^3 + 3fᵥ[1] * fᵥ[2] * e_f[2] + e_f[1] * fᵥ[3]
end
function (::DerivativeOfEnergy{4})(fᵥ, e_f)
    return e_f[4] * fᵥ[1]^4 +
           6fᵥ[1]^2 * fᵥ[2] * e_f[3] +
           (4fᵥ[1] * fᵥ[3] + 3fᵥ[3]^2) * e_f[2] +
           e_f[1] * fᵥ[4]
end

include("poly.jl")
