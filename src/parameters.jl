Base.@kwdef mutable struct Parameters
    Λ = 0.5879
    m0 = 0.005588
    Gs = 2.442/Λ^2 
    Gv = Gs/2
    eta_d = 0.75
    GD = eta_d*Gs
end

export Parameters
