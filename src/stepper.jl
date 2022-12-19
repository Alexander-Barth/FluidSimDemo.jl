function step(config,mask,uv,newuv,p)
    integrate!(config,mask,uv)
    p .= 0
    incompressibility!(config,mask,uv,p)
    boundary_conditions!(config,uv)
    advection!(config,mask,uv,newuv)
    boundary_conditions!(config,uv)
end


function main!(config,mask,uv,newuv,p,nmax; callback = nothing)
    for n = 1:nmax
        step(config,mask,uv,newuv,p)
        if !isnothing(callback)
            callback(mask,p,uv,n)
        end
    end
end
