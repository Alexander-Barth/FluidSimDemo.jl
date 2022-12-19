using FluidSimDemo
using BenchmarkTools

config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()

uv = (u,v)
newuv = (newu,newv)
nmax = 100

@btime FluidSimDemo.main!(config,mask,p,uv,newuv,nmax)

#@time FluidSimDemo.main!(config,mask,p,uv,newuv,nmax)
@show sum(p)
