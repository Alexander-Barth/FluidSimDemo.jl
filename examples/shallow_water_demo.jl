using FluidSimDemo
using Statistics
using BenchmarkTools


function sw_boundary_conditions!(config,(u,v))
    sz = (size(u,1)-1,size(u,2))
    @inbounds for j = 1:sz[2] # one single loop is marginally faster
        for i = 1:sz[1]
            if (i > 1)
                u[i,j] = u[i,j] * mask[i,j] * mask[i-1,j]
            end

            if (j > 1)
                v[i,j] = v[i,j] * mask[i,j] * mask[i,j-1]
            end
        end
    end
#    @show sz
end

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


x = config.Δxy[1] * (0:(size(mask,1)-1))
y = config.Δxy[2] * (0:(size(mask,2)-1))
L = 10 .* config.Δxy

courant_number = sqrt(config.grav*maximum(config.h)) * config.Δt / minimum(config.Δxy)

@show courant_number

p .= exp.(-(x .- mean(x)).^2/L[1]^2 .- ((y' .- mean(y)).^2/L[1]^2))
u .= 0
v .= 0

mask .= false
mask[2:end-1,2:end-1] .= true


using FluidSimDemo.Physics2D: interp, integrate!, incompressibility!, advection!#, free_surface!

#using FluidSimDemo.Physics2D: free_surface!
#using FluidSimDemo.PhysicsND: free_surface_nd!



#@btime free_surface!(config,mask,p,uv)

function step!(config,mask,p,uv,newuv)
    integrate!(config,mask,uv)
    config.boundary_conditions!(config,uv)
    FluidSimDemo.Physics2D.free_surface!(config,mask,p,uv)
    #FluidSimDemo.PhysicsND.free_surface!(config,mask,p,uv)
    config.boundary_conditions!(config,uv)
    advection!(config,mask,uv,newuv)
    config.boundary_conditions!(config,uv)
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
