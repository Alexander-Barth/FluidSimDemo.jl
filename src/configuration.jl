
"""
    set_mask!(config,mask,xy,(u,v))

Set mask based on the location of obstacle xy (tuple of floats).
Velocities on masked locations (land) are set to zero.
"""
function set_mask!(config,mask,xy,(u,v); radius = 0.15)
    sz = size(mask)
    h = config.h

    mask .= true
    @inbounds for j = 2:sz[2]-1
        for i = 2:sz[1]-1
            dx = (i-1 + 0.5) * h - xy[1]
            dy = (j-1 + 0.5) * h - xy[2]

            if dx^2 + dy^2 < radius^2
                mask[i,j] = false
                u[i,j] = 0
                u[i+1,j] = 0
                v[i,j] = 0
                v[i,j+1] = 0
            end
        end
    end

    mask[1,:] .= false
    mask[:,[1,end]] .= false
end




#=

+--------v[i,j+1]--------+
|                        |
|                        |
|                        |
u[i,j]      mask[i,j]     u[i+1,j]
|         div            |
|                        |
|                        |
+---------v[i,j]---------+

=#


function boundary_conditions!(config,(u,v))
    sz = (size(u,1)-1,size(u,2))

    # inflow
    @inbounds for j = 1:size(u,2)
        u[2,j] = config.u0
        # semi-lagrangian advection might query values on land
        u[1,j] = u[2,j]
    end

    @inbounds for i = 1:size(u,1)
        u[i,1] = u[i,2]
        u[i,end] = u[i,end-1]
    end

    @inbounds for j = 1:size(v,2)
        v[1,j] = v[2,j]
        v[end,j] = v[end-1,j]
    end
end


# sz size of the domain
# T type of variables
# xy initial position of obstace

function config_Karman_vortex_street(; sz = (300,100), T = Float32, xy = [0.4,0.5])
    # velocity components (staggered grid)
    u = zeros(T,sz[1]+1,sz[2])
    v = zeros(T,sz[1],sz[2]+1)

    # new velocity components
    newu = zeros(T,sz[1]+1,sz[2])
    newv = zeros(T,sz[1],sz[2]+1)

    # binary mask (1: sea, 0: land)
    mask = ones(Bool,sz)

    # pressure
    p = zeros(T,sz)

    config = (
        # inflow velocity
        u0 = T(2.), # m/s
        # grid resolution
        h = T(0.01), # m
        # time step
        Δt = T(1/60), # s
        # acceleration due to gravity
        g = (T(0.),T(0.)), # m/s²
        # density
        ρ = T(1000.), # kg/m³
        # overrelaxation parameter for Gauss-Seider
        overrelaxation = T(1.9),
        # iterations for the pressure solver
        iter_pressure = 40,
        # size of the domain
        sz = sz,
        xy = xy,
        boundary_conditions! = boundary_conditions!,
    )

    set_mask!(config,mask,xy,(u,v))
    boundary_conditions!(config,(u,v))

    return (config,mask,p,(u,v),(newu,newv))
end

