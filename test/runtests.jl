using Test

using FluidSimDemo
import FluidSimDemo.Physics2D
import FluidSimDemo.PhysicsND



config,mask,p,(u,v),(newu,newv) = FluidSimDemo.config_Karman_vortex_street()

F = [1 2; 3 4]
@test interp(F,(1,1)) ≈ 1.
@test interp(F,(1.5,1.5)) ≈ sum(F)/4
@test interp(F,(1.5,1)) ≈ (F[1,1] + F[2,1])/2
@test interp(F,(1.25,1)) ≈ (3*F[1,1] + F[2,1])/4
@test interp(F,(-1,1)) ≈ F[1,1]

# velocity components (staggered grid)

uv = (zeros(T,sz[1]+1,sz[2]),zeros(T,sz[1],sz[2]+1))
@time Physics2D.integrate!(config,mask,uv)

uvn = (zeros(T,sz[1]+1,sz[2]),zeros(T,sz[1],sz[2]+1))
@time PhysicsND.integrate!(config,mask,uvn)

@test uv[1] ≈ uvn[1]
@test uv[2] ≈ uvn[2]


#@btime Physics2D.integrate!(config,mask,uv)
#@btime PhysicsND.integrate4!(config,mask,uvn)


p = zeros(T,sz)
@time Physics2D.incompressibility!(config,mask,uv,p)

pn = zeros(T,sz)
@time PhysicsND.incompressibility!(config,mask,uvn,pn)

@test p ≈ pn



@time Physics2D.advection!(config,mask,uv,newuv);

newuvn = (zeros(T,sz[1]+1,sz[2]),zeros(T,sz[1],sz[2]+1))
@time PhysicsND.advection!(config,mask,uvn,newuvn);

@test uv[1] ≈ uvn[1]
@test uv[2] ≈ uvn[2]

@btime Physics2D.advection!(config,mask,uv,newuv);
@btime PhysicsND.advection!(config,mask,uv,newuv);




#p_save = copy(p)
#@show maximum(abs.(p - p_save))
#=

 gfortran  -O3 fluid_sim.f90 -o fluid_sim && ./fluid_sim
 mask                 28586
 cpu time   2.0293280000000000
 p   -466865.094
gher17:~/Julia/share (master)
$ gfortran  -O2 fluid_sim.f90 -o fluid_sim && ./fluid_sim
 mask                 28586
 cpu time   2.0689799999999998
 p   -466865.094



