function energy_q(q, param)
    return sqrt(q^2 + param.mq^2)
end

function energy_D(q, param)
    return sqrt(q^2 + param.mD^2)
end

function integrand_pol(T, μ, ω, param, p)
    return nothing 
end

function polarisation(T, μ, ω, param)
    return 4*param.m*integrate(p->p^2*integrand_pol(T,μ,ω,param,p),0, param.Λ)/(2π^2)
end

