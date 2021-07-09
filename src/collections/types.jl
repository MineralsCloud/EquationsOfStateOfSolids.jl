using AutoHashEquals: @auto_hash_equals
using ConstructionBase: constructorof
using EquationsOfState: Parameters, EquationOfState

export Murnaghan,
    Murnaghan1st,
    BirchMurnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    # PoirierTarantola4th,
    Vinet,
    EnergyEquation,
    PressureEquation,
    BulkModulusEquation

abstract type FiniteStrainParameters{N,T} <: Parameters{T} end
abstract type BirchMurnaghan{N,T} <: FiniteStrainParameters{N,T} end
abstract type PoirierTarantola{N,T} <: FiniteStrainParameters{N,T} end
abstract type Murnaghan{T} <: Parameters{T} end
"""
    Murnaghan1st(v0, b0, b′0, e0=zero(v0 * b0))

Create a Murnaghan first order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct Murnaghan1st{T} <: Murnaghan{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan1st{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    Murnaghan2nd(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Murnaghan second order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct Murnaghan2nd{T} <: Murnaghan{T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    Murnaghan2nd{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
"""
    BirchMurnaghan2nd(v0, b0, e0=zero(v0 * b0))

Create a Birch–Murnaghan second order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct BirchMurnaghan2nd{T} <: BirchMurnaghan{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
"""
    BirchMurnaghan3rd(v0, b0, b′0, e0=zero(v0 * b0))

Create a Birch–Murnaghan third order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.

!!! note
    The third-order equation (Equation (22)) becomes identical to the second-order equation
    when ``b′0 = 4`` (not ``0``!).
"""
@auto_hash_equals struct BirchMurnaghan3rd{T} <: BirchMurnaghan{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    BirchMurnaghan4th(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Birch–Murnaghan fourth order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.

!!! note
    The fourth-order equation becomes identical to the third-order equation when
    ```math
    b″0 = -\\frac{ 1 }{ 9b0 } (9b′0^2 - 63b′0 + 143).
    ```
"""
@auto_hash_equals struct BirchMurnaghan4th{T} <: BirchMurnaghan{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    BirchMurnaghan4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
"""
    PoirierTarantola2nd(v0, b0, e0=zero(v0 * b0))

Create a Poirier–Tarantola second order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola2nd{T} <: PoirierTarantola{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, e0)
end
"""
    PoirierTarantola3rd(v0, b0, b′0, e0=zero(v0 * b0))

Create a Poirier–Tarantola third order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola3rd{T} <: PoirierTarantola{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    PoirierTarantola4th(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Poirier–Tarantola fourth order equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct PoirierTarantola4th{T} <: PoirierTarantola{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    PoirierTarantola4th{T}(v0, b0, b′0, b″0, e0 = zero(v0 * b0)) where {T} =
        new(v0, b0, b′0, b″0, e0)
end
"""
    Vinet(v0, b0, b′0, e0=zero(v0 * b0))

Create a Vinet equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
@auto_hash_equals struct Vinet{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Vinet{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
@auto_hash_equals struct AntonSchmidt{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e∞::T
    AntonSchmidt{T}(v0, b0, b′0, e∞ = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e∞)
end
@auto_hash_equals struct Holzapfel{Z,T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    function Holzapfel{Z,T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {Z,T}
        @assert 1 <= Z <= 118 "elements are between 1 and 118!"
        return new(v0, b0, b′0, e0)
    end
end

function (::Type{T})(args...) where {T<:Parameters}
    E = Base.promote_typeof(args...)
    return constructorof(T){E}(args...)  # Cannot use `T.(args...)`! For `AbstractQuantity` they will fail!
end
function Murnaghan(args...)
    N = length(args)
    if N == 4
        return Murnaghan1st(args...)
    elseif N == 5
        return Murnaghan2nd(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end
"""
    BirchMurnaghan(args...)

Construct a `BirchMurnaghan` based on the length of arguments, where `e0` must be provided.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan3rd`](@ref), [`BirchMurnaghan4th`](@ref)
"""
function BirchMurnaghan(args...)
    N = length(args)
    if N == 3
        return BirchMurnaghan2nd(args...)
    elseif N == 4
        return BirchMurnaghan3rd(args...)
    elseif N == 5
        return BirchMurnaghan4th(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end
"""
    PoirierTarantola(args...)

Construct a `PoirierTarantola` based on the length of arguments, where `e0` must be provided.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref), [`PoirierTarantola4th`](@ref)
"""
function PoirierTarantola(args...)
    N = length(args)
    if N == 3
        return PoirierTarantola2nd(args...)
    elseif N == 4
        return PoirierTarantola3rd(args...)
    elseif N == 5
        return PoirierTarantola4th(args...)
    else
        throw(ArgumentError("unknown number of arguments $N."))
    end
end

"""
    EquationOfStateOfSolids{T<:Parameters}

Represent an equation of state of solids.
"""
abstract type EquationOfStateOfSolids{T} <: EquationOfState{T} end
"""
    EnergyEquation{T} <: EquationOfStateOfSolids{T}
    EnergyEquation(parameters::Parameters)

Construct an equation of state which evaluates the energy of the given `parameters`.
"""
struct EnergyEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end
EnergyEquation(eos::EquationOfStateOfSolids) = EnergyEquation(getparam(eos))
"""
    PressureEquation{T} <: EquationOfStateOfSolids{T}
    PressureEquation(parameters::Parameters)

Construct an equation of state which evaluates the pressure of the given `parameters`.
"""
struct PressureEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end
PressureEquation(eos::EquationOfStateOfSolids) = PressureEquation(getparam(eos))
"""
    BulkModulusEquation{T} <: EquationOfStateOfSolids{T}
    BulkModulusEquation(parameters::Parameters)

Construct an equation of state which evaluates the bulk modulus of the given `parameters`.
"""
struct BulkModulusEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end
BulkModulusEquation(eos::EquationOfStateOfSolids) = BulkModulusEquation(getparam(eos))
