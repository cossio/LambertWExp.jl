#= These functions give a fast but not very accurate 
approximation to W(e^x). I am not exporting these functions
because I am using the "exact" W(e^x) instead =#


"""fast approximation of W(e^x), with
max. abs. error 0.35 and max. rel. error 0.2.
It is concave, with continuous 2nd derivative
and discontinuous 3rd derivative at x=±2.
It is over 20x faster than LambertW.lambertw(exp(x))"""
function lambertwexp_approx(x::T) where {T<:Real}
    c0 = convert(float(T), 0.5895124241703698)
    c1 = convert(float(T), 0.2855227482747687)
    c2 = convert(float(T), 0.020207723988558527)
    c3 = convert(float(T), 0.0005797856921151388)
    c4 = convert(float(T), 0.00317192070085811)
    c5 = convert(float(T), 0.0003148433129769499)

    if x < -2
        exp(x)
    elseif x > 2
        x - log(x)
    else # -2 ≤ x ≤ 2 or NaN
        # interpolating polynomial with matching first and second derivatives
        #0.5895124241703698 + (0.2855227482747687 + (0.020207723988558527 + (0.0005797856921151388 + (0.00317192070085811 + 0.0003148433129769499x)x)x)x)x
        @evalpoly x c0 c1 c2 c3 c4 c5
    end
end


"""first derivative of the approximate lambertw(exp(x)) w.r.t x"""
function lambertwexp_approx_d1(x::T) where {T<:Real}
    c0 = convert(float(T), 0.2855227482747687)
    c1 = convert(float(T), 0.040415447977117054)
    c2 = convert(float(T), 0.0017393570763454165)
    c3 = convert(float(T), 0.01268768280343244)
    c4 = convert(float(T), 0.0015742165648847495)

    if x < -2
        exp(x)
    elseif x > 2
        one(x) - inv(x)
    else # -2 ≤ x ≤ 2 or NaN
        #0.2855227482747687 + (0.040415447977117054 + (0.0017393570763454165 + (0.01268768280343244 + 0.0015742165648847495x)x)x)x
        @evalpoly x c0 c1 c2 c3 c4
    end
end


"""second derivative of the approximate lambertw(exp(x)) w.r.t x"""
function lambertwexp_approx_d2(x::T) where {T<:Real}
    c0 = convert(float(T), 0.04041544797711706)
    c1 = convert(float(T), 0.0034787141526908326)
    c2 = convert(float(T), 0.03806304841029733)
    c3 = convert(float(T), 0.006296866259538999)

    if x < -2
        exp(x)
    elseif x > 2
        inv(x)^2
    else # -2 ≤ x ≤ 2 or NaN
        #0.04041544797711706 + (0.0034787141526908326 + (0.03806304841029733 + 0.006296866259538999x)x)x
        @evalpoly x c0 c1 c2 c3
    end
end
