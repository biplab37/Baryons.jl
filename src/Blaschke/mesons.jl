function imagpart_sigma(T, μ, ω, q, param)
    return nothing
end

function imagpart_sigma(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    if ω^2 < 4 * m^2
        return 0.0
    else
        return 3 * ω^2 * (1 - 4m^2 / ω^2) * sqrt(1 - 4m^2 / ω^2) *
               (1 - numberF(T, -μ, ω / 2) - numberF(T, μ, ω / 2)) / π
    end
end

function polarisation_sigma(T, μ, ω, q, param)
    return nothing
end

function polarisation_sigma2(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    integrand(Ep) = 12 * Ep * sqrt(Ep^2 - m^2) * (1 - m^2 / Ep^2) *
                    (
                        (1 - numberF(T, μ, Ep) - numberF(T, -μ, Ep)) *
                        (PrincipalValue(ω + 2Ep) - PrincipalValue(ω - 2Ep))
                    ) / π^2
    return integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function polarisation_sigma(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    cutoffE = sqrt(param.Λ^2 + m^2)
    integrand(Ep) = 12 * Ep * sqrt(Ep^2 - m^2) * (1 - m^2 / Ep^2) *
                    ((1 - numberF(T, μ, Ep) - numberF(T, -μ, Ep))) / (2π^2)
    return quadgk_cauchy(integrand, m, ω / 2, cutoffE) +
           quadgk_cauchy(integrand, m, -ω / 2, cutoffE)
end

function polarisation_sigma1(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    Ep(p) = sqrt(p^2 + m^2)
    integrand(p) = 12 * p^2 * (1 - m^2 / Ep(p)^2) *
                   (
                       (1 - numberF(T, μ, Ep(p)) - numberF(T, -μ, Ep(p))) *
                       (PrincipalValue(ω + 2Ep(p)) - PrincipalValue(ω - 2Ep(p)))
                   ) / π^2
    return integrate(integrand, 0, param.Λ)
end

function polarisation_sigma(ω, param)
    m = massgap(0.01, 0.0, param)[1]
    cutoffE = sqrt(param.Λ^2 + m^2)
    integrand(Ep) = 12 * Ep * sqrt(Ep^2 - m^2) * (1 - m^2 / Ep^2) / (2π^2)
    return quadgk_cauchy(integrand, m, ω / 2, cutoffE) +
           quadgk_cauchy(integrand, m, -ω / 2, cutoffE)
end

function mass_sigma_func(param)
    f(ω) = 1 / param.Gs - polarisation_sigma(ω, param)
    return f
end

function mass_sigma_func(T, μ, param)
    f(ω) = 1 / param.Gs - polarisation_sigma(T, μ, ω, param)
    return f
end

function mass_sigma(T, μ, param)
    f(ω) = 1 / param.Gs - polarisation_sigma(T, μ, ω, param)
    return fzero(f, 0.4)
end

function mass_sigma(trange::AbstractRange, μ, param)
    f(T, ω) = 1 / param.Gs - polarisation_sigma(T, μ, ω, param)
    masses = zeros(length(trange))

    for (i, t) in enumerate(trange)
        masses[i] = fzero(ω -> f(t, ω), 0.5)
    end
    return masses
end

function mass_sigma(T, μrange::AbstractRange, param)
    f(μ, ω) = 1 / param.Gs - polarisation_sigma(T, μ, ω, param)
    masses = zeros(length(μrange))

    for (i, μ) in enumerate(μrange)
        masses[i] = fzero(ω -> f(μ, ω), 0.5)
    end
    return masses
end

function phasesc_sigma(T, μ, ω, q, param)
    return phasesc(imagpart_sigma, T, μ, ω, q, param)
end

function imagpart_phi(T, μ, ω, q, param)
    return nothing
end

function imagpart_phi(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    if ω^2 < 4 * m^2
        return 0.0
    else
        return 3 * ω^2 * sqrt(1 - 4m^2 / ω^2) *
               (1 - numberF(T, -μ, ω / 2) - numberF(T, μ, ω / 2)) / π
    end
end

function polarisation_phi(T, μ, ω, q, param)
    return nothing
end

function polarisation_phi(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    cutoffE = sqrt(param.Λ^2 + m^2)
    integrand(Ep) =
        12 * Ep * sqrt(Ep^2 - m^2) * (1 - numberF(T, μ, Ep) - numberF(T, -μ, Ep)) / (2π^2)
    return quadgk_cauchy(integrand, m, ω / 2, cutoffE) +
           quadgk_cauchy(integrand, m, -ω / 2, cutoffE)
end

function mass_phi_func(T, μ, param)
    f(ω) = 1 / param.Gs - polarisation_phi(T, μ, ω, param)
    return f
end

function mass_phi(T, μ, param)
    f(ω) = 1 / param.Gs - polarisation_phi(T, μ, ω, param)
    return bisection(f, 0.0, 1.0)
end

function mass_phi(trange::AbstractRange, μ, param)
    f(T, ω) = 1 / param.Gs - polarisation_phi(T, μ, ω, param)
    masses = zeros(length(trange))

    for (i, t) in enumerate(trange)
        masses[i] = fzero(ω -> f(t, ω), 0.4)
    end
    return masses
end

function mass_phi(T, μrange::AbstractRange, param)
    f(μ, ω) = 1 / param.Gs - polarisation_phi(T, μ, ω, param)
    masses = zeros(length(μrange))

    for (i, μ) in enumerate(μrange)
        masses[i] = fzero(ω -> f(μ, ω), 0.4)
    end
    return masses
end

function phasesc_phi(T, μ, ω, q, param)
    return phasesc(imagpart_phi, T, μ, ω, q, param)
end

export mass_sigma, mass_phi, phasesc_sigma, phasesc_phi, mass_sigma_func, mass_phi_func
