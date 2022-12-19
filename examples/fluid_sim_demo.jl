using FluidSimDemo
using BenchmarkTools

config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()
myplot = FluidSimDemo.plotting(config,mask,p,(u,v), plot_every = 2)
nmax = 2000
FluidSimDemo.main!(config,mask,p,(u,v),(newu,newv),nmax,callback = myplot)


