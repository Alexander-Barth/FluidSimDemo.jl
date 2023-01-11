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

function incompressibility!(config,mask::AbstractArray{Bool,N},p,uv) where N
    Δt,ρ,h = config.Δt,config.ρ,config.h
    I = CartesianIndices(mask)
    Ifirst, Ilast = first(I), last(I)
    I1 = oneunit(Ifirst)
    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))

    cp = ρ * h / Δt

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
            for (u,uvec) in zip(uv,unitvecs)
                div += u[ij + uvec] - u[ij]
            end

            p_ = -div/nn
            p_ *= config.overrelaxation
            # pressure
            p[ij] += cp * p_

            for (u,uvec) in zip(uv,unitvecs)
                u[ij]      -= p_ *  mask[ij - uvec]
                u[ij+uvec] += p_ *  mask[ij + uvec]
            end
        end
    end
end

function advection!(config,mask::AbstractArray{Bool,N},uv::NTuple{N,AbstractArray{T,N}},newuv) where {N,T}
    Δt,h = config.Δt,config.h
    sz = size(mask)

    unitvecs = ntuple(i -> CartesianIndex(ntuple(==(i), Val(N))), Val(N))
    I = CartesianIndices(mask)
    Ifirst, Ilast = first(I), last(I)
    I1 = oneunit(Ifirst)

    ntuple(Val(N)) do i
        u = uv[i]
        newu = newuv[i]
        uvec = unitvecs[i]

        @inbounds for ij = Ifirst+I1:(Ilast-I1+uvec)
             if mask[ij] && mask[ij-unitvecs[i]]
                fij = ntuple(Val(N)) do j
                    @inbounds if i == j
                        ij[j] - (u[ij] * Δt/h)
                    else
                        v = uv[j]
                        vs = T(0.25) * (v[ij] + v[ij-uvec] +
                            v[ij+unitvecs[j]] + v[ij+unitvecs[j] - uvec])
                        ij[j] - (vs * Δt)/h
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
