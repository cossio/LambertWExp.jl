#= This code can be easily be extended to complex z
and other branches of W. I didn't do it because real
z was sufficient for my purposes. =#


export lambertwexp,
       lambertwexp_d1_from_W,
       lambertwexp_d2_from_W


"""Single update of the Fritsch iteration, w <- w + wε,
given the current iterate w and its residual z"""
function lambertw_fritsch_update(w::T, z::T) where T <: Real
    #= *******************************************************
    Fritsch, Shafer & Crowell 1973, https://doi.org/10.1145/361952.361970.
    See also Veberic 2012, https://doi.org/10.1016/j.cpc.2012.07.008,
    Sec. 2.3, which recommends this iteration for its fourth-order
    convergence (vs. third-order for Halley).

    For W(z) the iteration is w <- w(1 + ε), with

        zₙ = log(z / w) - w
        qₙ = 2(1 + w)(1 + w + (2/3)zₙ)
        ε  = (zₙ / (1 + w)) (qₙ - zₙ) / (qₙ - 2zₙ)

    ε is the relative correction applied to w, so it doubles as
    an error estimate.

    The update is computed as w + w * ε rather than w * (1 + ε):
    when ε is near roundoff, 1 + ε can round to 1 and lose the
    final ulp of correction that w + w * ε still delivers.
    ********************************************************** =#
    w1 = w + one(T)
    q = 2 * w1 * (w1 + z * convert(T, 2//3))
    ε = z / w1 * (q - z) / (q - 2 * z)
    return w + w * ε, abs(ε)
end


"""Fritsch iteration from initial guess w, where residual(w)
computes the residual zₙ of the current iterate"""
function lambertw_fritsch_iterate(residual, w::T, x::T, maxiter::Int) where T <: Real
    converged::Bool = false
    lastε = convert(T, Inf)
    for i in 1:maxiter
        wnew, ε = lambertw_fritsch_update(w, residual(w))
        #= wnew == w means w is a fixed point of the computed map;
        non-decreasing ε means the iteration reached its roundoff
        floor (e.g. a two-value cycle). Either way we are done. =#
        if wnew == w || ε ≥ lastε
            w = wnew
            converged = true
            break
        end
        w = wnew
        lastε = ε
    end
    converged || @warn "lambertwexp with z=", x, " did not converge in ", maxiter, " iterations."
    return w
end


"""Fritsch iteration for W(e^x), starting from an
initial guess appropriate to the magnitude of x"""
function lambertwexp_fritsch(x::T; maxiter::Int) where T <: Real
    if isnan(x)
        return x
    elseif !isfinite(x)
        return x > 0 ? x : zero(T)
    end

    if x > -2
        #= Residual zₙ = x - w - log(w), which never forms e^x and
        therefore cannot overflow. For x ≤ -2 this form is avoided:
        there log(w) ≈ x, so the roundoff of log(w) alone puts a
        floor of about eps * |x| on the achievable relative error. =#
        if x > 1
            # de Bruijn asymptotic expansion, W(e^x) ≈ x - log(x) + log(x)/x
            lx = log(x)
            w = x - lx + lx / x
        else
            w = lambertw_series_0(x)
        end
        return lambertw_fritsch_iterate(w -> x - w - log(w), w, x, maxiter)
    else
        #= Here a = e^x ≤ e⁻² is comfortably below overflow, so use the
        residual in its original form, zₙ = log(a / w) - w, which stays
        fully accurate as w -> 0 (a / w ≈ 1, no cancellation with x). =#
        a = exp(x)
        # W(a) = a - a² + (3/2)a³ - O(a⁴)
        w = a * (1 - a * (1 - convert(T, 3//2) * a))
        if a ≤ eps(one(T))
            #= The truncated series is already exact to roundoff.
            Returning here also avoids 0/0 in the residual when
            e^x underflows to zero. =#
            return w
        end
        return lambertw_fritsch_iterate(w -> log(a / w) - w, w, x, maxiter)
    end
end


"""Series expansion of W(e^z) about z = 0"""
function lambertw_series_0(x::T)::T where T<:Real
    0.5671432904097838 + (0.3618962566348892 +
                         (0.07367780517637275 +
                         (-0.001342859654990087 +
                         (-0.001636065147912496 + 0.00023214965556996332x)x)x)x)x
end


"""W(e^x), for real x and the principal branch of W"""
function lambertwexp(x::Real; maxiter::Integer = 1000)
    maxiter ≥ 0 || throw(ArgumentError("maxiter must be non-negative, got $maxiter"))
    return lambertwexp_fritsch(float(x); maxiter = Int(maxiter))
end


"""first derivative of W(e^x) with respect to x,
from the known value of W = W(e^x)"""
lambertwexp_d1_from_W(W::Real) = W / (1 + W)

"""second derivative of W(e^x) with respect to x,
from the known value of W = W(e^x)"""
lambertwexp_d2_from_W(W::Real) = W / (1 + W)^3
