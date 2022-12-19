function unstagger_vel((u,v),(u_r,v_r))
    sz = (size(u,1)-1,size(u,2))

    @inbounds for j = 1:sz[2]
        for i = 1:sz[1]
            u_r[i,j] = 0.5 * (u[i,j] + u[i+1,j])
            v_r[i,j] = 0.5 * (v[i,j] + v[i,j+1])
        end
    end
end



