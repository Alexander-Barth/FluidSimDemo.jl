"""
Model based on on the Navier-Stokes equations:

    ∂ₜ u + u ∇u = -∇p/ρ + g
          ∇ ⋅ u = 0

where u is the n-dimensional velocity (typically 2 or 3), ρ density,
g the accerlation vector due to gravity, p the pressure and ∇ is the nabla
operator. The density ρ is assumed constant.
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
