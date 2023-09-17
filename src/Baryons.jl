module Baryons

using UsefulFunctions
using QuadGK

include("helperfunctions.jl")
include("parameters.jl")
include("Blaschke/gapeqn.jl")
include("Blaschke/polarisation.jl")
include("Wang/gapeqn.jl")
include("Wang/polarisation.jl")
include("Wang/spectralfunction.jl")

end
