using FluidSimDemo

# setup an the example confguration of the Karman vortex street
config,mask,p,(u,v) = FluidSimDemo.config_Karman_vortex_street()
# p is the pressure
# u,v is the x- and y- velocity

# create a plotting call-back function
myplot = FluidSimDemo.plotting(config,mask,p,(u,v), plot_every = 2)
# number of time steps
nmax = 2000
FluidSimDemo.main!(config,mask,p,(u,v),nmax,callback = myplot)


