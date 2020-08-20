module Fitting

using ConstructionBase: constructorof, setproperties
using LsqFit: curve_fit, coef
using PolynomialRoots: roots
using Polynomials: Polynomial, fit, derivative, coeffs, derivative
using Serialization: serialize
using Unitful: AbstractQuantity, ustrip, unit

using ..Collections:
    EquationOfStateOfSolids,
    FiniteStrainParameters,
    Parameters,
    PressureEOS,
    EnergyEOS,
    BulkModulusEOS,
    orderof,
    strain_from_volume,
    volume_from_strain,
    whatstrain

export linfit, nonlinfit

_islocalminimum(y, x) = derivative(y, 2)(x) > 0  # If 2nd derivative at `x > 0`, `x` is a local minimum.

function _findlocalminima(y)
    y′ = derivative(y, 1)
    pool = real(filter(isreal, roots(coeffs(y′))))  # Complex volumes are meaningless
    return [x for x in pool if _islocalminimum(y, x)]
end

_findminimum(y) = _findminimum(y, _findlocalminima(y))
function _findminimum(y, localminima)  # Find the minimal in the minima
    # https://stackoverflow.com/a/21367608/3260253
    if isempty(localminima)
        error("no real local minima found!")  # For some polynomials, could be all complex
    else
        y0, i = findmin(y.(localminima))
        x0 = localminima[i]
        return x0, y0
    end
end

function linfit(eos::EnergyEOS{<:FiniteStrainParameters}, volumes, energies)
    deg = orderof(eos)
    v0_init = iszero(eos.param.v0) ? volumes[findmin(energies)[2]] : eos.param.v0
    st = whatstrain(eos.param)
    strains = map(strain_from_volume(st, v0_init), volumes)
    poly = fit(strains, energies, deg)
    f0, e0 = _findminimum(poly)
    v0_final = volume_from_strain(st, v0_init)(f0)
    f′ᵥ = -(2f0 + 1)^(5 / 2) / 3 / v0_final
    f″ᵥ = (2f0 + 1)^4 / (5 / (9 * v0_final^2))
    f‴ᵥ = -(2f0 + 1)^(11 / 2) * (40 / (27 * v0_final^3))
    e″ᵥ = derivative(poly, 2)(f0) * f′ᵥ^2 + derivative(poly, 1)(f0) * f″ᵥ
    b0 = v0_final * e″ᵥ
    e‴ᵥ =
        f′ᵥ^3 * derivative(poly, 3) +
        3f′ᵥ * f″ᵥ * derivative(poly, 2) +
        derivative(poly, 1) * f‴ᵥ
    b′0 = -v0_final * e‴ᵥ / e″ᵥ - 1
    return constructorof(eos.param)(v0_final, b0, b′0, e0)
end
    e″ᵥ = derivative(poly, 2)(f0) * f′ᵥ^2 + derivative(poly, 1)(f0) * f″ᵥ
    b0 = v0_final * e″ᵥ
    e‴ᵥ =
        f′ᵥ^3 * derivative(poly, 3) +
        3f′ᵥ * f″ᵥ * derivative(poly, 2) +
        derivative(poly, 1) * f‴ᵥ
    b′0 = -v0_final * e‴ᵥ / e″ᵥ - 1
    return constructorof(eos.param)(v0_final, b0, b′0, e0)
end

_energy′ᵥ(fᵥ, e_f) = e_f[1] * fᵥ[1]
_energy″ᵥ(fᵥ, e_f) = e_f[2] * fᵥ[1]^2 + e_f[1] * fᵥ[2]
_energy‴ᵥ(fᵥ, e_f) = e_f[3] * fᵥ[1]^3 + 3fᵥ[1] * fᵥ[2] * e_f[2] + e_f[1] * fᵥ[3]
_energy⁗ᵥ(fᵥ, e_f) =
    e_f[4] * fᵥ[1]^4 +
    6fᵥ[1]^2 * fᵥ[2] * e_f[3] +
    (4fᵥ[1] * fᵥ[3] + 3fᵥ[3]^2) * e_f[2] +
    e_f[1] * fᵥ[4]

function nonlinfit(
    eos::EquationOfStateOfSolids,
    xs,
    ys;
    xtol = 1e-8,
    gtol = 1e-2,
    maxiter::Integer = 1000,
    min_step_quality = 1e-3,
    good_step_quality = 0.75,
    silent = true,
    saveto = "",
)
    model = createmodel(eos)
    p0, xs, ys = _prepare(eos, xs, ys)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        xs,
        ys,
        collect(last.(p0));
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = !silent,
    )
    result = if fit.converged
        constructorof(typeof(eos.param))(
            map(coef(fit), first.(p0), _mapfields(unit, eos.param)) do x, c, u
                x / c * u
            end,
        )
    else
        @error "fitting is not converged, change initial parameters!"
        if !isinteractive() && isempty(saveto)
            saveto = string(rand(UInt)) * ".jls"
        end
        nothing
    end
    if !isempty(saveto)
        _savefit(saveto, fit)
    end
    _checkresult(result)
    return result
end

function createmodel(::S) where {T,S<:EquationOfStateOfSolids{T}}  # Do not export!
    constructor = constructorof(S) ∘ constructorof(T)
    return (x, p) -> map(constructor(p), x)
end

function _savefit(file, fit)  # Do not export!
    open(file, "w") do io
        @info "saving raw fitted data to '$file'..."
        serialize(io, fit)
    end
end

function _checkresult(param::Parameters)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "fitted `v0` or `b0` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

function _prepare(eos, xs, ys)  # Do not export!
    xs, ys = _collect_float(xs), _collect_float(ys)  # `xs` & `ys` may not be arrays
    if eos isa EnergyEOS && iszero(eos.param.e0)
        eos = EnergyEOS(setproperties(eos.param; e0 = minimum(ys)))  # Energy minimum as e0
    end
    return _ustrip_all(eos, xs, ys)
end

# No need to constrain `eltype`, `ustrip` will error if `Real` and `AbstractQuantity` are met.
function _ustrip_all(eos, xs, ys)  # Do not export!
    xs, ys = ustrip.(unit(eos.param.v0), xs), ustrip.(_yunit(eos), ys)
    punit = unit(eos.param.e0) / unit(eos.param.v0)
    return map(fieldnames(typeof(eos.param))) do f
        x = getfield(eos.param, f)
        if f == :b0
            ustrip(punit, oneunit(x)) => ustrip(punit, float(x))
        elseif f == :b″0
            ustrip(punit^(-1), oneunit(x)) => ustrip(punit^(-1), float(x))
        elseif f == :b‴0
            ustrip(punit^(-2), oneunit(x)) => ustrip(punit^(-2), float(x))
        else
            1 => ustrip(float(x))
        end
    end, xs, ys
end
_yunit(eos::EnergyEOS) = unit(eos.param.e0)
_yunit(eos::Union{PressureEOS,BulkModulusEOS}) = unit(eos.param.e0) / unit(eos.param.v0)

_collect_float(x) = collect(float.(x))  # Do not export!

_mapfields(f, x) = (f(getfield(x, i)) for i in 1:nfields(x))  # Do not export!

Base.float(p::Parameters) = constructorof(typeof(p))(_mapfields(float, p)...)  # Not used here but may be useful

end
