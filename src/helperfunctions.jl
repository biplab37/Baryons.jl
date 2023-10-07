function integrate(func, a, b)
    return quadgk(func, a, b, rtol = 1e-3, maxevals=1e5)[1]
end

function En(p,m)
    return sqrt(p^2 + m^2)
end

function fzero(f, guess)
    return nlsolve(x->f(x...), [guess]).zero[1]
end

function cauchy_quadgk(g, a, b; kws...)
    a < 0 < b || throw(ArgumentError("domain must include 0"))
    g₀ = g(0)
    g₀int = b == -a ? zero(g₀) : g₀ * log(abs(b/a)) / (b - a)
    return quadgk_count(x -> (g(x)-g₀)/x + g₀int, a, 0, b; kws...)
end
