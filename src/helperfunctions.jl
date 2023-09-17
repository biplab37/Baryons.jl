function integrate(func, a, b)
    return quadgk(func, a, b, maxiter = 10000)[1]
end
