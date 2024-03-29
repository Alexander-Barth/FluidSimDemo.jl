function step!(config,mask,p,uv,newuv)
    integrate!(config,mask,uv)
    config.boundary_conditions!(config,mask,uv)
    p .= 0
    incompressibility!(config,mask,p,uv)
    config.boundary_conditions!(config,mask,uv)
    advection!(config,mask,uv,newuv)
    config.boundary_conditions!(config,mask,uv)
end


function main!(config,mask,p,uv,newuv,nmax::Integer;
               callback = nothing,
               )
    for n = 1:nmax
        #@show uv[1][1:30,50]
        step!(config,mask,p,uv,newuv)
        if !isnothing(callback)
            callback(mask,p,uv,n)
        end
    end
end

function main!(config,mask,p,uv,nmax::Integer; callback = nothing,
               newuv = similar.(uv),
               boundary_conditions! = nothing,
               )

    for u in newuv
        u .= 0
    end

    main!(config,mask,p,uv,newuv,nmax; callback = callback)
end
