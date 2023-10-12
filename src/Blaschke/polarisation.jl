function energy_q(T, μ, q, param)
    return sqrt(q^2 + massgap(T, μ, param)[1]^2)
end

function energy_D(T, μ, q, param)
    return sqrt(q^2 + mass_diquark(T, μ, param)^2)
end

function integrand_pol(T, μ, ω, param, p)
    eq = energy_q(T, μ, p, param)
    eD = energy_D(T, μ, p, param)
    term1 = (1 - numberF(T, μ, eq) + numberB(T, μ, eD)) / (ω + eq + eD)
    term2 = (1 - numberF(T, -μ, eq) + numberB(T, -μ, eD)) * PrincipalValue(ω - eq - eD)
    term3 = (-numberF(T, -μ, eq) - numberB(T, μ, eD)) * PrincipalValue(ω - eq + eD)
    term4 = (numberF(T, μ, eq) + numberB(T, -μ, eD)) * PrincipalValue(ω + eq - eD)
    return (term1 - term2 + term3 + term4) / (4 * eq * eD)
end

function polarisation_baryon(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    return 4 * m * integrate(p -> p^2 * integrand_pol(T, μ, ω, param, p), 0, param.Λ) /
           (2π^2)
end
