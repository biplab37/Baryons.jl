function integrate(func, a, b)
    return quadgk(func, a, b, rtol=1e-3, maxevals=1e5)[1]
end

function En(p, m)
    return sqrt(p^2 + m^2)
end

function fzero(f, guess)
    return nlsolve(x -> f(x...), [guess]).zero[1]
end

function quadgk_cauchy(f, a, c, b)
    if (a - c) * (b - c) >= 0
        return integrate(x -> f(x) / (x - c), a, b)
    else
        fc = f(c)
        g(x) = (f(x) - fc) / (x - c)
        return quadgk(g, a, c, b)[1] + fc * log(abs((b - c) / (a - c)))
    end
end

function deriv(f, point)
    return central_fdm(5, 1)(f, point)
end
