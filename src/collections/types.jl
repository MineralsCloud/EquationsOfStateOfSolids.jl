abstract type EquationOfStateOfSolidsParameters{T} end
const Parameters = EquationOfStateOfSolidsParameters
abstract type FiniteStrainParameters{N,T} <: Parameters{T} end
abstract type BirchMurnaghan{N,T} <: FiniteStrainParameters{N,T} end
abstract type PoirierTarantola{N,T} <: FiniteStrainParameters{N,T} end
"""
    Murnaghan(v0, b0, b′0, e0=zero(v0 * b0))

Create a Murnaghan equation of state.

This equation of state can have units. The units are specified in
[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s `@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
"""
abstract type Murnaghan{T} <: Parameters{T} end
@auto_hash_equals struct Murnaghan1st{T} <: Murnaghan{T}
    v0::T
    b0::T
    b′0::T
    e0::T
    Murnaghan1st{T}(v0, b0, b′0, e0 = zero(v0 * b0)) where {T} = new(v0, b0, b′0, e0)
end
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
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `e0`: the energy of solid at zero pressure.
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
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure. "`″`" can be typed by `\\pprime<tab>`.
- `e0`: the energy of solid at zero pressure.
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
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
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
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
- `b″0`: the second-order pressure-derivative bulk modulus of solid at zero pressure. "`″`" can be typed by `\\pprime<tab>`.
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
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure. "`′`" can be typed by `\\prime<tab>`.
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

abstract type EquationOfState{T<:Parameters} end
abstract type EquationOfStateOfSolids{T} <: EquationOfState{T} end
struct EnergyEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct PressureEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end
struct BulkModulusEquation{T} <: EquationOfStateOfSolids{T}
    param::T
end

abstract type EquationFrom{T<:EquationOfStateOfSolids} end
struct EnergyFrom{T} <: EquationOfStateOfSolids{T}
    eos::T
end
struct PressureFrom{T} <: EquationOfStateOfSolids{T}
    eos::T
end
struct BulkModulusFrom{T} <: EquationOfStateOfSolids{T}
    eos::T
end
(x::EnergyFrom)(v) = EnergyEquation(getparam(x.eos))(v)
(x::EnergyFrom{<:EnergyEquation})(v) = x.eos(v)
(x::PressureFrom)(v) = PressureEquation(getparam(x.eos))(v)
(x::PressureFrom{<:PressureEquation})(v) = x.eos(v)
(x::BulkModulusFrom)(v) = BulkModulusEquation(getparam(x.eos))(v)
(x::BulkModulusFrom{<:BulkModulusEquation})(v) = x.eos(v)
