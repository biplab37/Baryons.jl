Base.@kwdef mutable struct Parameters
    Λ = 0.5879
    m0 = 0.005588
    Gs = 2.442/Λ^2 
    Gv = 1.221/Λ^2
end

export Parameters
