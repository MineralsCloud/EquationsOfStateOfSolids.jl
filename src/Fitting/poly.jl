using Polynomials: derivative, coeffs
using PolynomialRoots: roots

struct GlobalMinimumNotFound
    msg::String
end

# `x` is a real number and `y` is a real polynomial.
function islocalminimum(x, y)
    y″ = derivative(y, 2)
    y″ₓ = y″(x)
    return y″ₓ > zero(y″ₓ)  # If 2nd derivative at `x` > 0, `(x, y(x))` is a local minimum.
end

# Return the real local minima of a real polynomial
function localminima(y, root_thr=1e-20)  # `y` is a polynomial (could be complex).
    y′ = derivative(y, 1)
    allroots = roots(coeffs(y′); polish=true, epsilon=root_thr)
    realroots = real.(filter(isreal, allroots))  # Complex volumes are meaningless.
    if isempty(realroots)
        throw(GlobalMinimumNotFound("consider changing `root_thr`?"))
    end
    localminimaₓ = filter(Base.Fix2(islocalminimum, y), realroots)
    if isempty(localminimaₓ)
        throw(GlobalMinimumNotFound("all local extrema are local maxima!"))
    end
    return zip(localminimaₓ, map(y, localminimaₓ))
end

# https://stackoverflow.com/a/21367608/3260253
function globalminimum(y, root_thr=1e-20)  # Find the minimum of the local minima
    @assert isreal(y)
    localminimaₓ = first(localminima(y, root_thr))
    y₀, index = findmin(last(localminima(y, root_thr)))
    x₀ = localminimaₓ[index]
    return x₀, y₀
end
