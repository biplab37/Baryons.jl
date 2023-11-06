module Baryons

using UsefulFunctions
# using FiniteDifferences, ForwardDiff
using QuadGK, NLsolve

include("helperfunctions.jl")
include("parameters.jl")
include("polarisation_common.jl")

include("Blaschke/gapeqn.jl")
include("Blaschke/mesons.jl")
include("Blaschke/diquark.jl")
include("Blaschke/polarisation.jl")
include("Blaschke/baryon_mass.jl")

include("Wang/gapeqn.jl")
include("Wang/polarisation.jl")
include("Wang/spectralfunction.jl")

end
