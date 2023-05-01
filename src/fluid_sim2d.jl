module Physics2D


#=

    +--------v[i,j+1]--------+
    |                        |
    |                        |
    |         p[i,j]         |
 u[i,j]      mask[i,j]     u[i+1,j]
    |          div           |
    |                        |
    |                        |
    +---------v[i,j]---------+

=#

function interp(F,(fi,fj))
    fi = clamp(fi,1,size(F,1))
    fj = clamp(fj,1,size(F,2))

    i = min(floor(Int,fi),size(F,1)-1)
    j = min(floor(Int,fj),size(F,2)-1)

    a = fi - i
    b = fj - j

    return @inbounds (((1-a)*((1-b) * F[i,j]   + b * F[i,j+1]) +
                           a*((1-b) * F[i+1,j] + b * F[i+1,j+1])))
end

function integrate!(config,mask,(u,v))
    Δt,g = config.Δt,config.g
    sz = size(mask)

    @inbounds for j = 1:sz[2] # one single loop is marginally faster
        for i = 1:sz[1]
            if (i > 1)
                #u[i,j] += g[1] * Δt * mask[i,j] * mask[i-1,j]
                u[i,j] = (u[i,j] + g[1] * Δt) * mask[i,j] * mask[i-1,j]
            end

            if (j > 1)
                #v[i,j] += g[2] * Δt * mask[i,j] * mask[i,j-1]
                v[i,j] = (v[i,j] + g[2] * Δt) * mask[i,j] * mask[i,j-1]
            end
        end
    end
end

function free_surface!(config,mask,η,(u,v))
    h = config.h[1,1] # depth m
    g = config.grav
    Δt = config.Δt
    Δx, Δy = config.Δxy
    inv_Δx = 1/Δx
    inv_Δy = 1/Δy

    imax,jmax = size(η)
    Δt_over_Δx = Δt/Δx
    Δt_over_Δy = Δt/Δy
    g_Δt_over_Δx = g*Δt/Δx
    g_Δt_over_Δy = g*Δt/Δy

    @inbounds for j = 1:jmax-1
        for i = 1:imax-1
            η[i,j] = η[i,j] - (  (h*u[i+1,j] - h*u[i,j]) * Δt_over_Δx
                               + (h*v[i,j+1] - h*v[i,j]) * Δt_over_Δy)
        end
    end

    # single loop is a bit faster
    @inbounds for j = 1:jmax
        for i = 1:imax
            if i > 2
                u[i,j] = u[i,j] - (η[i,j] - η[i-1,j]) * g_Δt_over_Δx
            end

            if j > 2
                v[i,j] = v[i,j] - (η[i,j] - η[i,j-1]) * g_Δt_over_Δy
            end
        end
    end
end

function incompressibility!(config,mask,p,(u,v))
    Δt,ρ = config.Δt,config.ρ
    Δx, Δy = config.Δxy
    sz = size(mask)

    inv_Δx = 1/Δx
    inv_Δy = 1/Δy
    ΔA = Δx * Δy
    cp = ρ / Δt

    # Gauss-Seidel
    @inbounds for iter = 1:config.iter_pressure
        for j = 2:sz[2]-1
            for i = 2:sz[1]-1

                if !mask[i,j]
                    continue
                end

                #  number of direct neightbors with water
                nn = mask[i+1,j] + mask[i-1,j] + mask[i,j+1] + mask[i,j-1]

                if nn == 0
                    continue
                end

                div = inv_Δx * (u[i+1,j] - u[i,j]) + inv_Δy * (v[i,j+1] - v[i,j])
                p_ = -(div * ΔA)/nn
                p_ *= config.overrelaxation
                # pressure
                p[i,j] += cp * p_

                u[i,j]   -= p_ * mask[i-1,j] * inv_Δx
                u[i+1,j] += p_ * mask[i+1,j] * inv_Δx
                v[i,j]   -= p_ * mask[i,j-1] * inv_Δy
                v[i,j+1] += p_ * mask[i,j+1] * inv_Δy
            end
        end
    end
end

function advection!(config,mask,(u,v),(newu,newv))
    T = eltype(u)
    Δt = config.Δt
    Δx, Δy = config.Δxy
    sz = size(mask)
    inv_Δx = 1/Δx
    inv_Δy = 1/Δy

    @inbounds for j = 2:sz[2]-1
        for i = 2:sz[1]
            if mask[i,j] && mask[i-1,j]
                v_u = T(0.25) * (v[i,j] + v[i,j+1] + v[i-1,j] + v[i-1,j+1])
                fi = i + (-u[i,j] * Δt) * inv_Δx
                fj = j + (-v_u * Δt) * inv_Δy
                newu[i,j] = interp(u,(fi,fj))
            end
        end
    end

    @inbounds for j = 2:sz[2]
        for i = 2:sz[1]-1
            if mask[i,j] && mask[i,j-1]
                u_v = T(0.25) * (u[i,j] + u[i+1,j] + u[i,j-1] + u[i+1,j-1])
                fi = i + (-u_v * Δt) * inv_Δx
                fj = j + (-v[i,j] * Δt) * inv_Δy
                newv[i,j] = interp(v,(fi,fj))
            end
        end
    end

    u .= newu
    v .= newv
end

end
