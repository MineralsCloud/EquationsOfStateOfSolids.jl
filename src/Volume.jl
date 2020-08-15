module Volume

using UnPack: @unpack

using ..Collections

export volumeof

function volumeof(eos::PressureEoss{<:Murnaghan}, p)
    @unpack v0, b0, b′0, e0 = eos.param
    return v0 * (1 + b′0 / b0 * p)^(-1 / b′0)
end
function volumeof(eos::EnergyEoss{<:BirchMurnaghan3rd}, e)
    @unpack v0, b0, b′0, e0 = eos.param
    fs = [
        -0.3333333333333333 * b0 * v0 * Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) +
        Complex(-1.8898815748423097, 3.273370907915164) *
        Power(b0, 2) *
        Power(v0, 2) *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            -1 // 3,
        ) +
        Complex(-0.01469815788859444, -0.02545795624071486) *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            1 // 3,
        ),
        -0.3333333333333333 * b0 * v0 * Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) +
        Complex(-1.8898815748423097, -3.273370907915164) *
        Power(b0, 2) *
        Power(v0, 2) *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            -1 // 3,
        ) +
        Complex(-0.01469815788859444, 0.02545795624071486) *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            1 // 3,
        ),
        -0.3333333333333333 * b0 * v0 * Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) +
        3.7797631496846193 *
        Power(b0, 2) *
        Power(v0, 2) *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            -1 // 3,
        ) +
        0.02939631577718888 *
        Power(-4.0 * b0 * v0 + b0 * b′0 * v0, -1) *
        Power(
            69984.0 * Power(b0, 2) * e * Power(v0, 2) +
            -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
            4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
            -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
            34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
            -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
            -1458.0 * Power(b0, 3) * Power(v0, 3) +
            Power(
                -2.125764e6 * Power(b0, 6) * Power(v0, 6) + Power(
                    69984.0 * Power(b0, 2) * e * Power(v0, 2) +
                    -34992.0 * Power(b0, 2) * b′0 * e * Power(v0, 2) +
                    4374.0 * Power(b0, 2) * Power(b′0, 2) * e * Power(v0, 2) +
                    -69984.0 * Power(b0, 2) * e0 * Power(v0, 2) +
                    34992.0 * Power(b0, 2) * b′0 * e0 * Power(v0, 2) +
                    -4374.0 * Power(b0, 2) * Power(b′0, 2) * e0 * Power(v0, 2) +
                    -1458.0 * Power(b0, 3) * Power(v0, 3),
                    2,
                ),
                1 // 2,
            ),
            1 // 3,
        ),
    ]
    v = map(volume_from_strain(Eulerian(), v0) ∘ Complex, fs)
    return v
end

function Power(f, g)
    typeof(f) == Complex{Float64} ? fc = f : fc = Complex(f, 0)
    if imag(fc) == 0.0
        fc = real(fc) + 0.0im
    end
    fc^g
end

end
