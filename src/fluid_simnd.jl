module PhysicsND

function interp(F::AbstractArray{T,N},fij) where {T,N}
    fij_ = clamp.(fij,1,size(F))
    ij = min(floor.(Int,fij_),size(F) .- 1)
    a = fij_ .- ij

    Finterp = zero(T)
    @inbounds for ind in CartesianIndices(ntuple(i -> 0:1,Val(N)))
        w = one(T)
        for side in 1:N
            if ind[side] == 0
                w *= (1-a[side])
            else
                w *= a[side]
            end
        end
        Finterp += w * F[ind + CartesianIndex(ij)]
    end
    return Finterp
end

function integrate!(config,mask::AbstractArray{Bool,N},uv) where N
    Δt,g = config.Δt,config.g
    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))

    @inbounds for ij in CartesianIndices(mask)
        ntuple(Val(N)) do i
            @inbounds if ij[i] > 1
                #uv[i][ij] += g[i] * Δt * mask[ij] * mask[ij-unitvecs[i]]
                uv[i][ij] = (uv[i][ij] + g[i] * Δt) * mask[ij] * mask[ij-unitvecs[i]]
            end
        end
    end
end

function free_surface!(config,mask::AbstractArray{Bool,N},η,uv) where N
    #T = eltype(uv[1])
    h = config.h # depth m
    h_uv = config.h_uv # depth m
    grav = config.grav
    Δt = config.Δt
    Δxy = config.Δxy
    sz = size(mask)

    Δt_over_Δxy = Δt ./ Δxy
    g_Δt_over_Δxy = grav * Δt ./ Δxy

    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))
    I = CartesianIndices(mask)
    Ifirst, Ilast = first(I), last(I)
    I1 = oneunit(Ifirst)

    # compute surface elevation based on diverge of the transport
    @inbounds for ij = Ifirst:Ilast-I1
        for (u,h_u,Δt_over_Δx,uvec) in zip(uv,h_uv,Δt_over_Δxy,unitvecs)
            η[ij] -= (h_u[ij + uvec]*u[ij + uvec] - h_u[ij]*u[ij]) * Δt_over_Δx
        end
        # same speed
        # ntuple(Val(N)) do i
        #     @inbounds begin
        #         u = uv[i]
        #         h_u = h_uv[i]
        #         Δt_over_Δx = Δt_over_Δxy[i]
        #         uvec = unitvecs[i]
        #         η[ij] -= (h_u[ij + uvec]*u[ij + uvec] - h_u[ij]*u[ij]) * Δt_over_Δx
        #     end
        # end

    end

    # update velocities based on pressure gradient
    # single loop is a bit faster
    @inbounds for ij in I
        #significantly slower
        # for (i,(u,g_Δt_over_Δx,uvec)) in enumerate(zip(uv,g_Δt_over_Δxy,unitvecs))
        #     if ij[i] > 2
        #         u[ij] -= (η[ij] - η[ij - uvec]) * g_Δt_over_Δx
        #     end
        # end
        ntuple(Val(N)) do i
            @inbounds begin
                u = uv[i]
                if ij[i] > 2
                    u[ij] -= (η[ij] - η[ij - unitvecs[i]]) * g_Δt_over_Δxy[i]
                end
            end
        end
    end
end

function incompressibility!(config,mask::AbstractArray{Bool,N},p,uv) where N
    Δt,ρ = config.Δt,config.ρ
    Δxy = config.Δxy

    I = CartesianIndices(mask)
    Ifirst, Ilast = first(I), last(I)
    I1 = oneunit(Ifirst)
    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))

    inv_Δxy = 1 ./ Δxy
    ΔA = prod(Δxy)
    cp = ρ / Δt

    # Gauss-Seidel
    @inbounds for iter = 1:config.iter_pressure
        for ij = Ifirst+I1:Ilast-I1
            if !mask[ij]
                continue
            end

            #  number of direct neightbors with water
            nn = 0
            for uvec in unitvecs
                nn += mask[ij+uvec]+mask[ij-uvec]
            end

            if nn == 0
                continue
            end

            div = zero(eltype(uv[1]))
            for (u,uvec,inv_Δx) in zip(uv,unitvecs,inv_Δxy)
                div += (u[ij + uvec] - u[ij]) * inv_Δx
            end

            p_ = -(div * ΔA)/nn
            p_ *= config.overrelaxation
            # pressure
            p[ij] += cp * p_

            for (u,uvec,inv_Δx) in zip(uv,unitvecs,inv_Δxy)
                u[ij]      -= p_ *  mask[ij - uvec] * inv_Δx
                u[ij+uvec] += p_ *  mask[ij + uvec] * inv_Δx
            end
        end
    end
end

function advection!(config,mask::AbstractArray{Bool,N},uv::NTuple{N,AbstractArray{T,N}},newuv) where {N,T}
    Δt = config.Δt
    Δxy = config.Δxy
    sz = size(mask)

    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))
    I = CartesianIndices(mask)
    Ifirst, Ilast = first(I), last(I)
    I1 = oneunit(Ifirst)

    inv_Δxy = 1 ./ Δxy

    ntuple(Val(N)) do i
        u = uv[i]
        newu = newuv[i]
        uvec = unitvecs[i]

        @inbounds for ij = Ifirst+I1:(Ilast-I1+uvec)
             if mask[ij] && mask[ij-unitvecs[i]]
                fij = ntuple(Val(N)) do j
                    @inbounds if i == j
                        ij[j] - (u[ij] * Δt) * inv_Δxy[j]
                    else
                        v = uv[j]
                        vs = T(0.25) * (v[ij] + v[ij-uvec] +
                            v[ij+unitvecs[j]] + v[ij+unitvecs[j] - uvec])
                        ij[j] - (vs * Δt) * inv_Δxy[j]
                    end
                end
                 newu[ij] = interp(u,fij)
             end
        end
    end

    ntuple(Val(N)) do i
        u = uv[i]
        newu = newuv[i]
        @inbounds for ij in eachindex(u)
            u[ij] = newu[ij]
        end
    end
end
end
