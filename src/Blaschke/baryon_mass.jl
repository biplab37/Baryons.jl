function coupling(T, μ, param)
    return 4 * g2(T, μ, param) / massgap(T, μ, param)[1]
end

function g2(T, μ, param)
    mD = mass_diquark(T, μ, param)
    return 4 * mD / deriv(x -> polarisation_diquark(T, μ, x, param), mD)
end

function mass_baryon(T, μ, param; guess=0.9)
    return fzero(x -> 1 - coupling(T, μ, param) * polarisation_baryon(T, μ, x, param), guess)
end

function mass_baryon_func(T, μ, param)
    return x -> 1 - coupling(T, μ, param) * polarisation_baryon(T, μ, x, param)
end

function mass_baryon(trange::AbstractRange, μ, param; initial_guess=1.0)
    f(T, ω) = 1 - coupling(T, μ, param) * polarisation_baryon(T, μ, ω, param)
    masses = zeros(length(trange))
    guess = initial_guess
    for (i, t) in enumerate(trange)
        guess = fzero(ω -> f(t, ω), guess)
        masses[i] = guess
    end
    return masses
end

function mass_baryon(T, μrange::AbstractRange, param; initial_guess=0.9)
    f(μ, ω) = 1 - coupling(T, μ, param) * polarisation_baryon(T, μ, ω, param)
    masses = zeros(length(μrange))
    guess = initial_guess
    for (i, μ) in enumerate(μrange)
        guess = fzero(ω -> f(μ, ω), guess)
        masses[i] = guess
    end
    return masses
end

export mass_baryon
