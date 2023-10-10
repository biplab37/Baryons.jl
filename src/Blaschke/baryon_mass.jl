function coupling(param)
    return 4 * g2(param) / param.m
end

function g2(param)
    return 4 * param.mD / derivative(polarisation, param.mD)
end

function baryon_mass(T, μ, param)
    return bisection(x -> 1 - coupling(param) * polarisation(T, μ, x, param), 0, 1.2)
end
