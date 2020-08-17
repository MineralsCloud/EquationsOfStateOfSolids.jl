module Fitting

using ConstructionBase: constructorof, setproperties
using LsqFit: curve_fit, coef
using Serialization: serialize
using Unitful: AbstractQuantity, ustrip, unit

using ..Collections:
    EquationOfStateOfSolids,
    FiniteStrainParameters,
    Parameters,
    PressureEOS,
    EnergyEOS,
    BulkModulusEOS

export linfit, nonlinfit

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
    p0, xs, ys, rules = _preprocess(eos, xs, ys)
    fit = curve_fit(  # See https://github.com/JuliaNLSolvers/LsqFit.jl/blob/f687631/src/levenberg_marquardt.jl#L3-L28
        model,
        xs,
        ys,
        p0;
        x_tol = xtol,
        g_tol = gtol,
        maxIter = maxiter,
        min_step_quality = min_step_quality,
        good_step_quality = good_step_quality,
        show_trace = !silent,
    )
    result = if fit.converged
        constructorof(typeof(eos))(map(coef(fit), rules) do x, (from, to)
            x * to |> from
        end)
    else
        if !isinteractive()
            saveto = string(rand(UInt)) * ".jls"
        end
        @error "fitting is not converged, change initial parameters!"
    end
    if !isempty(saveto)
        _savefit(saveto, fit)
    end
    _checkparam(result)
    return result
end # function nonlinfit

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

function _checkparam(param::FiniteStrainParameters)  # Do not export!
    if param.v0 <= zero(param.v0) || param.b0 <= zero(param.b0)
        @error "fitted `v0` or `b0` is negative!"
    end
    # if PressureEoss(param)(minimum(v)) >= param.b0
    #     @warn "use higher order EOS!"
    # end
end

_collect_float(x) = collect(float.(x))  # Do not export!

function _preprocess(eos, xs, ys)  # Do not export!
    xs, ys = _collect_float(xs), _collect_float(ys)  # `xs` & `ys` may not be arrays
    if eos isa EnergyEOS && iszero(eos.param.e0)
        eos = EnergyEOS(setproperties(eos.param; e0 = minimum(ys)))  # Energy minimum as e0
    end
    return _ustrip_all(eos, xs, ys)
end

# No need to constrain `eltype`, `ustrip` will error if `Real` and `AbstractQuantity` are met.
function _ustrip_all(eos::EnergyEOS, vs, es)  # Do not export!
    vunit, eunit = unit(eos.param.v0), unit(eos.param.e0)
    punit = eunit / vunit
    vs, es = ustrip.(vunit, vs), ustrip.(eunit, es)
    rules = map(fieldnames(typeof(eos.param))) do f
        from = unit(getfield(eos.param, f))
        to = if f == :b0
            @set! eos.param.b0 = ustrip(punit, eos.param.b0)
            punit
        elseif f == :b′0
            1
        elseif f == :b′′0
            @set! eos.param.b′′0 = ustrip(punit^(-1), eos.param.b′′0)
            punit^(-1)
        elseif f == :b′′′0
            @set! eos.param.b′′′0 = ustrip(punit^(-2), eos.param.b′′′0)
            punit^(-2)
        elseif f == :v0
            vunit
        elseif f == :e0
            eunit
        end
        from => to
    end
    return collect(_splat(ustrip ∘ float, eos.param)), vs, es, rules
end
function _ustrip_all(eos::Union{PressureEOS,BulkModulusEOS}, vs, ps)
    vunit, eunit = unit(eos.param.v0), unit(eos.param.e0)
    punit = eunit / vunit
    vs, ps = ustrip.(vunit, vs), ustrip.(punit, ps)
    rules = map(fieldnames(typeof(eos.param))) do f
        from = unit(getfield(eos.param, f))
        to = if f == :b0
            @set! eos.param.b0 = ustrip(punit, eos.param.b0)
            punit
        elseif f == :b′0
            1
        elseif f == :b′′0
            @set! eos.param.b′′0 = ustrip(punit^(-1), eos.param.b′′0)
            punit^(-1)
        elseif f == :b′′′0
            @set! eos.param.b′′′0 = ustrip(punit^(-2), eos.param.b′′′0)
            punit^(-2)
        elseif f == :v0
            vunit
        elseif f == :e0
            eunit
        end
        from => to
    end
    return collect(_splat(ustrip ∘ float, eos.param)), vs, ps, rules
end

_mapfields(f, x) = (f(getfield(x, i)) for i in 1:nfields(x))  # Do not export!

Base.float(p::Parameters) = constructorof(typeof(p))(_mapfields(float, p)...)  # Not used here but may be useful

end
