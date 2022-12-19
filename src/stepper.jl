function step(config,mask,p,uv,newuv)
    integrate!(config,mask,uv)
    p .= 0
    incompressibility!(config,mask,p,uv)
    boundary_conditions!(config,uv)
    advection!(config,mask,uv,newuv)
    boundary_conditions!(config,uv)
end


function main!(config,mask,p,uv,newuv,nmax; callback = nothing)
    for n = 1:nmax
        step(config,mask,p,uv,newuv)
        if !isnothing(callback)
            callback(mask,p,uv,n)
        end
    end
end
