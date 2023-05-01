function unstagger_vel((u,v),(u_r,v_r))
    sz = (size(u,1)-1,size(u,2))

    @inbounds for j = 1:sz[2]
        for i = 1:sz[1]
            u_r[i,j] = 0.5 * (u[i,j] + u[i+1,j])
            v_r[i,j] = 0.5 * (v[i,j] + v[i,j+1])
        end
    end
end


function stagger_r2uv!(h,(h_u,h_v))
    sz = size(h)

    @inbounds for j = 1:sz[2]+1
        for i = 1:sz[1]+1
            if j <= sz[2]
                h_u[i,j] = (h[min(i,sz[1]),j] + h[max(i-1,1),j]) / 2
            end

            if i <= sz[1]
                h_v[i,j] = (h[i,min(j,sz[2])] + h[i,max(j-1,1)]) / 2
            end
        end
    end

    return (h_u,h_v)
end


function stagger_r2uv(h)
    sz = size(h)
    h_u = similar(h,(sz[1]+1,sz[2]))
    h_v = similar(h,(sz[1],sz[2]+1))
    return stagger_r2uv!(h,(h_u,h_v))
end
