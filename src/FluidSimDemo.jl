module FluidSimDemo

include("fluid_sim_makie.jl")
include("fluid_simnd.jl")
include("fluid_sim2d.jl")


#using .PhysicsND: interp, integrate!, incompressibility!, advection!
using .Physics2D: interp, integrate!, incompressibility!, advection!

include("configuration.jl")
include("stepper.jl")

include("utils.jl")

end # module FluidSimDemo
