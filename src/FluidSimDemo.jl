"""
Model based on on the Navier-Stokes equations:

    ∂ₜ u + u ∇u = -∇p/ρ
          ∇ ⋅ u = 0

The density ρ is assumed constant.
"""
module FluidSimDemo

include("fluid_simnd.jl")
include("fluid_sim2d.jl")


using .PhysicsND: interp, integrate!, incompressibility!, advection!
#using .Physics2D: interp, integrate!, incompressibility!, advection!

include("configuration.jl")
include("stepper.jl")

include("utils.jl")

end # module FluidSimDemo
