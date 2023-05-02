using FluidSimDemo
using Statistics
using BenchmarkTools

using FluidSimDemo:sw_boundary_conditions!, sw_initial_conditions!

# setup an the example confguration of the Karman vortex street
config,mask,p,(u,v) = FluidSimDemo.config_Karman_vortex_street(
    Δt = 2e-3,u0 = 0.)
# p is the pressure
# u,v is the x- and y- velocity

sz = config.sz
T = eltype(config.Δt)

h = ones(T,config.sz)
h_uv = FluidSimDemo.stagger_r2uv(h)

config = (config..., boundary_conditions! = sw_boundary_conditions!,
          grav = T(9.81),
          h = h,
          h_uv = h_uv,
          )

# create a plotting call-back function
#myplot = FluidSimDemo.plotting(config,mask,p,(u,v), plot_every = 2)
myplot = nothing

# number of time steps
#nmax = 2000
nmax = 2000
nmax = 20



courant_number = sqrt(config.grav*maximum(config.h)) * config.Δt / minimum(config.Δxy)

sw_initial_conditions!(config,mask,p,(u,v))

@show courant_number



using FluidSimDemo.Physics2D: interp, integrate!, incompressibility!, advection!#, free_surface!

#using FluidSimDemo.Physics2D: free_surface!
#using FluidSimDemo.PhysicsND: free_surface_nd!



#@btime free_surface!(config,mask,p,uv)

function step!(config,mask,p,uv,newuv)
    #integrate!(config,mask,uv)
    config.boundary_conditions!(config,mask,uv)
    FluidSimDemo.Physics2D.free_surface!(config,mask,p,uv)
    config.boundary_conditions!(config,mask,uv)
    advection!(config,mask,uv,newuv)
    config.boundary_conditions!(config,mask,uv)
end

uv = (u,v)
newuv = similar.(uv)
callback = myplot

for n = 1:nmax
    #@show uv[1][1:30,50]
    step!(config,mask,p,uv,newuv)
    if !isnothing(callback)
        callback(mask,p,uv,n)
    end
end


using Test
@test u[150,50] ≈ -0.090978034f0
@test p[150,50] ≈ -0.17512816f0
