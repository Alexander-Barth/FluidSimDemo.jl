using FluidSimDemo
using BenchmarkTools

config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()

T = eltype(u)
uv = (u,v)
newuv = (newu,newv)
nmax = 1

# warmup
FluidSimDemo.main!(config,mask,p,uv,newuv,1)

config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()
uv = (u,v)
newuv = (newu,newv)
nmax = 100

@time FluidSimDemo.main!(config,mask,p,uv,newuv,nmax)

@show sum(p)
