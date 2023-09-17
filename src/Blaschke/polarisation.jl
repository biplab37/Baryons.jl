function energy_q(T,\mu,q, param)
    return sqrt(q^2 + mass_q(T, μ, param)^2)
end

function energy_D(T,μ,q, param)
    return sqrt(q^2 + mass_D(T,μ,param)^2)
end

function integrand_pol(T, μ, ω, param, p)
    eq = energy_q(q, param)
    eD = energy_D(q, param)
    term1 = (1 - numberF(T, μ, q) + numberB(T, μ, eD))/(ω + eq + eD)
    term2 = (1 - numberF(T, -μ, q) + numberB(T, -μ, eD))*PrincipalValue(ω - eq - eD)
    term3 = (-numberF(T, -μ, q) - numberB(T, μ, eD))*PrincipalValue(ω - eq + eD)
    term4 = (numberF(T, μ, q) + numberB(T, -μ, eD))*PrincipalValue(ω + eq - eD)
    return (term1 + term2 + term3 + term4)/(4*eq*eD)
end

function polarisation(T, μ, ω, param)
    return 4*param.m*integrate(p->p^2*integrand_pol(T,μ,ω,param,p),0, param.Λ)/(2π^2)
end

