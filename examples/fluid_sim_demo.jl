using FluidSimDemo
using BenchmarkTools



config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()



FluidSimDemo.boundary_conditions!(config,(u,v))

myplot = FluidSimDemo.fluidplot_setup(config,mask,p,(u,v))


#s = fluidplot_setup(cx,cy,mask,p,(u,v),0)


uv = (u,v)
newuv = (newu,newv)
nmax = 10

# 0.0135 s
#@time main!(config,mask,uv,newuv,p,1)
#FluidSimDemo.main!(config,mask,uv,newuv,p,100)


@time FluidSimDemo.main!(config,mask,uv,newuv,p,100, callback =
    (mask,p,uv,n) -> begin
                         if n % 2 == 0
                             myplot(mask,p,uv,n)
                         end
                         end
    
    )

#@time main!(config,mask,uv,newuv,p,nmax)



@show sum(p)
