function energy_q(T, μ, q, param)
    return sqrt(q^2 + massgap(T, μ, param)[1]^2)
end

function energy_D(T, μ, q, param)
    return sqrt(q^2 + mass_diquark(T, μ, param)^2)
end

#TODO: rewrite this function to handle the principal value integral better.
function integrand_pol(T, μ, ω, param, p)
    eq = energy_q(T, μ, p, param)
    eD = energy_D(T, μ, p, param)
    term1 = (1 - numberF(T, μ, eq) + numberB(T, μ, eD)) / (ω + eq + eD)
    term2 = (1 - numberF(T, -μ, eq) + numberB(T, -μ, eD)) * PrincipalValue(ω - eq - eD)
    term3 = (-numberF(T, -μ, eq) - numberB(T, μ, eD)) * PrincipalValue(ω - eq + eD)
    term4 = (numberF(T, μ, eq) + numberB(T, -μ, eD)) * PrincipalValue(ω + eq - eD)
    return (term1 - term2 + term3 + term4) / (4 * eq * eD)
end

function polarisation_baryon1(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    return 4 * m * integrate(p -> p^2 * integrand_pol(T, μ, ω, param, p), 0, param.Λ) / (2π^2)
end

function integrand_baryon(T, μ, ω, param, p, mq, mD, pole_p)
    eq = sqrt(p^2 + mq^2)
    eD = sqrt(p^2 + mD^2)
    term1 = (1 - numberF(T, μ, eq) + numberB(T, μ, eD)) * (p - pole_p) / (ω + eq + eD)
    term2 = (1 - numberF(T, -μ, eq) + numberB(T, -μ, eD)) * (p - pole_p) / (ω - eq - eD)
    term3 = (-numberF(T, -μ, eq) - numberB(T, μ, eD)) * (p - pole_p) / (ω - eq + eD)
    term4 = (numberF(T, μ, eq) + numberB(T, -μ, eD)) * (p - pole_p) / (ω + eq - eD)
    return (term1 - term2 + term3 + term4) / (4 * eq * eD)
end

function two_body_pole(ω, mq, mD)
    deter = (1 - (mq + mD)^2 / ω^2) * (1 - (mq - mD)^2 / ω^2)
    if deter < 0
        return 0.0
    else
        return ω * sqrt(deter) / 2
    end
end

function polarisation_baryon(T, μ, ω, param)
    mq = massgap(T, μ, param)[1]
    mD = mass_diquark(T, μ, param)
    pole_p = two_body_pole(ω, mq, mD)

    return 4 * mq * quadgk_cauchy(p -> p^2 * integrand_baryon(T, μ, ω, param, p, mq, mD, pole_p), 0, pole_p, param.Λ) / (2π^2)
end