using LambertW


export lambertwexp

#= This can easily be extended to complex z and
other branches of W. =#


"""Halley's iteration for W(e^z), with initial guess x0"""
function lambertwexp_halley(z::T, x0::T; maxiter::Int) where T <: Number
    #= *******************************************************
    Based on https://github.com/jlapeyre/LambertW.jl.git,
    which computes W(x). We here needed to compute W(e^x).

    Halley's iteration for W(z) is given in:

    Corless, et al 1996, https://doi.org/10.1007/BF02124750

    In Eq. 5.9 of that paper, substitute z by e^z, and
    divide numerator and denominator of the fraction by
    e^z. This gives the iteration used here for W(e^z).

    This routine is best when z ≥ 0. If z < 0 it might be better
    to call lambert(exp(z)) directly.
    ********************************************************** =#

    two_t = convert(T, 2)
    x = x0
    lastx = x
    lastdiff = zero(T)
    converged::Bool = false
    for i in 1:maxiter
        ex = exp(x - z)
        xexz = x * ex - 1
        x1 = x + 1
        x -= xexz / (ex * x1 - (x + two_t) * xexz / (two_t * x1 ))
        xdiff = abs(lastx - x)
        if xdiff <= 3 * eps(abs(lastx)) || lastdiff == xdiff  # second condition catches two-value cycle
            converged = true
            break
        end
        lastx = x
        lastdiff = xdiff
    end
    converged || @warn "lambertwexp with z=", z, " did not converge in ", maxiter, " iterations."
    return x
end


"""Provides initial guesses for Halley's iteration 
to compute W(e^z), using the principal branch of W"""
function lambertw_branch_zero_exp(x::T; maxiter::Integer)::T where T<:Real
    if isnan(x)
        return x
    elseif !isfinite(x)
        return x > 0 ? x : zero(x)
    # elseif iszero(x)
    #     # omega constant, solution of xe^x=1
    #     return T(0.567143290409783872999968662210355)
    end

    one_t = one(T)
    itwo_t = 1 / convert(T, 2)

    if x > zero(T)
        lx = x
        llx = log(lx)
        x0 = lx - llx - log(one_t - llx / lx) * itwo_t
    else
        x0 = (567//1000) * exp(x)
    end

    return lambertwexp_halley(x, x0; maxiter = maxiter)
end


"""W(e^x), for real x and the principal branch of W"""
function lambertwexp(x::Real; maxiter::Integer = 1000)
    maxiter ≥ 0 || throw(ArgumentError("maxiter must be non-negative, got $maxiter"))
    if x < 0
        lambertw(exp(x))
    else
        lambertw_branch_zero_exp(float(x); maxiter = maxiter)
    end
end
