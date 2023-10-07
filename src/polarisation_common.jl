## real part from krammers-kronig relation

function realpart(imagpart::Function, T, μ, ω, q, param::Parameters)
    integrand(v) = 2*v* imagpart(T, μ, ω, q, param)* (PrincipalValue(v^2 - ω^2) - PrincipalValue(v^2))/π
    integrand_inv(v) = integrand(1/(1 - v))/(1 - v)^2
    return integrate(integrand, 0., 1.) + integrate(integrand_inv, 0., 1.)
end

function realpart(imagpart::Function, T, μ, ω, param::Parameters)
    integrand(v) = 2*v* imagpart(T, μ, ω, param)* (PrincipalValue(v^2 - ω^2) - PrincipalValue(v^2))/π
    integrand_inv(v) = integrand(1/(1 - v))/(1 - v)^2
    return integrate(integrand, 0., 1.) + integrate(integrand_inv, 0., 1.)
end

function phasesc(imagpart::Function, T, μ, ω, q, param::Parameters)
    impi = imagpart(T, μ, ω, q, param)
    repi = realpart(imagpart, T, μ, ω, q, param)
    return angle(Complx(repi, -impi))
end

function mass(polarisation, T, μ, coupling, param::Parameters; initial_guess = 0.5)
    spectral(ω) = 1/coupling - polarisation(T, μ, ω, param)
    return bisection(spectral, 0.0, 1.0)
end
