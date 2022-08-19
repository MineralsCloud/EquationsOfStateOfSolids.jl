using ConstructionBase: constructorof
using EquationsOfState: Parameters, EquationOfState
using Functors: @functor
using StructHelpers: @batteries

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

The energy and pressure equations are:
```math
\\begin{align*}
    E(V) &= E_0+B_0 V_0\\left[\\frac{1}{B_0'\\left(B_0'-1\\right)}\\left(\\frac{V}{V_0}\\right)^{1-B_0'}+\\frac{V}{B_0' V_0}-\\frac{1}{B_0'-1}\\right]\\\\
    P(V) &= \\frac{B_0}{B_0'}\\left[\\left(\\frac{V_0}{V}\\right)^{B_0'}-1\\right]
\\end{align*}
```

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
struct Murnaghan1st{T} <: Murnaghan{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan1st{T}(v0, b0, b′0, e0=zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
"""
    Murnaghan2nd(v0, b0, b′0, b″0, e0=zero(v0 * b0))

Create a Murnaghan second order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
struct Murnaghan2nd{T} <: Murnaghan{T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    function Murnaghan2nd{T}(v0, b0, b′0, b″0, e0=zero(v0 * b0)) where {T}
        return new(v0, b0, b′0, b″0, e0)
    end
end
"""
    BirchMurnaghan2nd(v0, b0, e0=zero(v0 * b0))

Create a Birch–Murnaghan second order equation of state.

The energy, pressure, and bulk modulus equations are:
```math
\\begin{align*}
    E(V) &= E_0 + \\frac{9}{8} B_0 V_0 \\left(\\left( V / V_0 \\right)^{-2 / 3}-1\\right)^{2}\\\\
    P(V) &= \\frac{3}{2} B_{0}\\left(x^{-7 / 3}-x^{-5 / 3}\\right)\\\\
    B(V) &= B_0(7f+1)(2f+1)^{5 / 2}
\\end{align*}
```
where ``x = V / V_0``, and ``f = \\frac{ 1 }{ 2 } \\bigg[ \\Big( \\frac{ V_0 }{ V } \\Big)^{2/3} - 1 \\bigg]``.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure.
"""
struct BirchMurnaghan2nd{T} <: BirchMurnaghan{2,T}
    v0::T
    b0::T
    e0::T
    BirchMurnaghan2nd{T}(v0, b0, e0=zero(v0 * b0)) where {T} = new(v0, b0, e0)
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
struct BirchMurnaghan3rd{T} <: BirchMurnaghan{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    BirchMurnaghan3rd{T}(v0, b0, b′0, e0=zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
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
    B''_0 = -\\frac{ 1 }{ 9B_0 } (9B'_0^2 - 63B'_0 + 143).
    ```
"""
struct BirchMurnaghan4th{T} <: BirchMurnaghan{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    function BirchMurnaghan4th{T}(v0, b0, b′0, b″0, e0=zero(v0 * b0)) where {T}
        return new(v0, b0, b′0, b″0, e0)
    end
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
struct PoirierTarantola2nd{T} <: PoirierTarantola{2,T}
    v0::T
    b0::T
    e0::T
    PoirierTarantola2nd{T}(v0, b0, e0=zero(v0 * b0)) where {T} = new(v0, b0, e0)
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
struct PoirierTarantola3rd{T} <: PoirierTarantola{3,T}
    v0::T
    b0::T
    b′0::T
    e0::T
    PoirierTarantola3rd{T}(v0, b0, b′0, e0=zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
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
struct PoirierTarantola4th{T} <: PoirierTarantola{4,T}
    v0::T
    b0::T
    b′0::T
    b″0::T
    e0::T
    function PoirierTarantola4th{T}(v0, b0, b′0, b″0, e0=zero(v0 * b0)) where {T}
        return new(v0, b0, b′0, b″0, e0)
    end
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
struct Vinet{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Vinet{T}(v0, b0, b′0, e0=zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
struct AntonSchmidt{T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e∞::T
    AntonSchmidt{T}(v0, b0, b′0, e∞=zero(v0 * b0)) where {T} = new(v0, b0, b′0, e∞)
end
struct Holzapfel{Z,T} <: Parameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    function Holzapfel{Z,T}(v0, b0, b′0, e0=zero(v0 * b0)) where {Z,T}
        @assert 1 <= Z <= 118 "elements are between 1 and 118!"
        return new(v0, b0, b′0, e0)
    end
end

@batteries Murnaghan1st eq = true hash = true
@batteries Murnaghan2nd eq = true hash = true
@batteries BirchMurnaghan2nd eq = true hash = true
@batteries BirchMurnaghan3rd eq = true hash = true
@batteries BirchMurnaghan4th eq = true hash = true
@batteries PoirierTarantola2nd eq = true hash = true
@batteries PoirierTarantola3rd eq = true hash = true
@batteries PoirierTarantola4th eq = true hash = true
@batteries Vinet eq = true hash = true
@batteries AntonSchmidt eq = true hash = true
@batteries Holzapfel eq = true hash = true
@functor Murnaghan1st
@functor Murnaghan2nd
@functor BirchMurnaghan2nd
@functor BirchMurnaghan3rd
@functor BirchMurnaghan4th
@functor PoirierTarantola2nd
@functor PoirierTarantola3rd
@functor PoirierTarantola4th
@functor Vinet
@functor AntonSchmidt
@functor Holzapfel

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
