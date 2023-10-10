function imagpart_diquark(T, μ, ω, q, param)
    return nothing
end

function imagpart_diquark(T, μ, ω, param)
    return nothing
end

function polarisation_diquark(T, μ, ω, q, param)
    return nothing
end

function polarisation_diquark1(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    integrand(Ep) = -(16 / π^2) * Ep * sqrt(Ep^2 - m^2) * ((Ep - μ) * (1 - 2 * numberF(T, μ, Ep)) * PrincipalValue(ω^2 - 4 * (Ep - μ)^2) + (Ep + μ) * (1 - 2 * numberF(T, -μ, Ep)) * PrincipalValue(ω^2 - 4 * (Ep + μ)^2))
    return integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function polarisation_diquark(T, μ, ω, param)
    m = massgap(T, μ, param)[1]
    integrand(Ep) = -(16 / π^2) * Ep * sqrt(Ep^2 - m^2) * ((Ep - μ) * (1 - 2 * numberF(T, μ, Ep)) * PrincipalValue(ω^2 - 4 * (Ep - μ)^2) + (Ep + μ) * (1 - 2 * numberF(T, -μ, Ep)) * PrincipalValue(ω^2 - 4 * (Ep + μ)^2))
    return integrate(integrand, m, sqrt(param.Λ^2 + m^2))
end

function mass_diquark_func(T, μ, param)
    f(ω) = 1 / param.GD - polarisation_diquark(T, μ, ω, param)
    return f
end

function mass_diquark(T, μ, param)
    f(ω) = 1 / param.GD - polarisation_diquark(T, μ, ω, param)
    return fzero(f, 0.5)
end

function mass_diquark(trange::AbstractRange, μ, param; initial_guess=0.5)
    f(T, ω) = 1 / param.GD - polarisation_diquark(T, μ, ω, param)
    masses = zeros(length(trange))
    guess = initial_guess
    for (i, t) in enumerate(trange)
        guess = fzero(ω -> f(t, ω), guess)
        masses[i] = guess
    end
    return masses
end

function mass_diquark(T, μrange::AbstractRange, param)
    f(μ, ω) = 1 / param.GD - polarisation_diquark(T, μ, ω, param)
    masses = zeros(length(μrange))
    guess = initial_guess
    for (i, μ) in enumerate(μrange)
        guess = fzero(ω -> f(μ, ω), 0.0, 1.0)
        masses[i] = guess
    end
    return masses
end

function phasesc_diquark(T, μ, ω, q, param)
    return phasesc(imagpart_diquark, T, μ, ω, q, param)
end

export mass_diquark, mass_diquark_func
