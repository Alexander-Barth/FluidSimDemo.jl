
This package solves the inviscit, incompressible Navier Stokes equations in 2 or 3 (or any other) dimensions.



https://user-images.githubusercontent.com/9881475/212127012-f02edcfb-7b9c-4c51-9056-cac4a3718876.mp4



# Installation

In a julia terminal, run the following commands to install this package:

```julia
using Pkg
Pkg.add(url="https://github.com/Alexander-Barth/FluidSimDemo.jl")
```

Then run the [examples/fluid_sim_demo.jl](examples/fluid_sim_demo.jl).


# Mathematical background

The evolution equation of the velocity $\mathbf u$ is given by:

$$
\frac{d \mathbf u}{dt} + \mathbf u \cdot \nabla \mathbf u = -\frac{1}{\rho} \nabla p + \mathbf g
$$

where $p$ is the pressure, $\rho$ is the constant density and $\mathbf g$ and external force like gravity. $\nabla$ is the [nabla operator](https://en.wikipedia.org/wiki/Del).

The flow is subjected to the following incompressibility constraint:

$$
\nabla \cdot \mathbf u = 0
$$


These equations are solved in this different steps:

1. Apply external forces $\mathbf g$ (in particular gravity)

$$
\mathbf {u'}^{(n)} = \mathbf u^{(n-1)} + \Delta t \mathbf g
$$

Impose on boundaries $\mathbf {u'}^{(n)} = 0$ or the corresponding inflow/outflow velocity. 


2. Advection of velocity field  $\mathbf {u''}^{(n)}$. The combined force of inertia and external forces acting from $n-1$ to $n$ are written as $\mathbf F$:


$$
\mathbf {u''}^{(n)} = \mathbf {u}^{(n-1)} + \Delta t \mathbf F
$$

3. The pressure is computed by requiring that the divergence of velocity must remain zero:

$$
\nabla \cdot \left( \mathbf F - \frac{1}{\rho} \nabla p^{(n)}  \right) = 0
$$

which leads to:

$$
\nabla^2 p^{(n)} = \rho  \nabla \cdot \mathbf {F}^{(n)}
$$

Note that it is not necessary to compute $\mathbf F$ separatly as its divergence can be compute from 
$\mathbf {u''}^{(n)}$.


$$
\nabla \cdot \mathbf {F}^{(n)} = \frac{1}{\Delta t} \nabla \cdot \mathbf {u''}^{(n)}
$$


$$
\nabla \cdot \mathbf {u''}^{(n)} = \Delta t  \nabla \cdot \mathbf F
$$

Choose pessure such that:


$$
\nabla \cdot \left( \mathbf {u''}^{(n)}  - \Delta t \frac{1}{\rho} \nabla p^{(n)}  \right) = 0
$$


$$
\nabla^2 p^{(n)} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf {u''}^{(n)}
$$


The pressure is solved iteratively (using Gauss-Seidel with overrelaxation) with a fixed number of iterations. In the following algorithm $\leftarrow$ is the assignement operator. For simplicity, the algorithm is outlined for the 2D case:
Intitialize the iteration:

$$
\begin{alignat*}{2}
& u'''\_{i,j} &&\leftarrow u''^{(n)}\_{i,j} \\
& v'''\_{i,j} &&\leftarrow v''^{(n)}\_{i,j} \\
& p\_{i,j}    &&\leftarrow 0 \\
\end{alignat*}
$$

The time index $n$ is dropped as all parmeters are from the same time instance.

Compute pressure adjustement by:

$$
\Delta p \leftarrow -\frac{\rho \Delta x \Delta y}{4 \Delta t} \left(\frac{u'''\_{i+1,j} - u'''\_{i,j}}{\Delta x} + \frac{v'''\_{i,j+1} - v'''\_{i,j}}{\Delta y} \right)
$$

Adjust the pressure

$$
p_{i,j} \leftarrow p_{i,j} + \Delta p
$$

Update the velocity accordingly

$$\begin{array}{cc}
u'''\_{i,j} &\leftarrow u'''\_{i,j} - \frac{\Delta p\_{i,j} - \Delta p\_{i-1,j}}{\Delta x} \frac{\Delta t}{\rho} \\
v'''\_{i,j} &\leftarrow v'''\_{i,j} - \frac{\Delta p\_{i,j} - \Delta p\_{i,j-1}}{\Delta y} \frac{\Delta t}{\rho} \\
\end{array}
$$

After convergence, the velocity $u'''\_{i,j}$ and $v'''\_{i,j}$ correspond to the velocity for the next time step $\mathbf {u}^{(n)}$. On can show 
that $\mathbf {u}^{(n)}$ and $p$ satisfy the following equations:

$$
\mathbf {u}^{(n)} = \mathbf {u}^{(n-1)} + \Delta t \mathbf F - \frac{\Delta t}{\rho} \nabla p
$$
