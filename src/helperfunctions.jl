function integrate(func, a, b)
    return quadgk(func, a, b, rtol = 1e-2)[1]
end

function En(p,m)
    return sqrt(p^2 + m^2)
end
