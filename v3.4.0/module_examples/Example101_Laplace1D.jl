#=

# 101: 1D Laplace equation 
([source code](@__SOURCE_URL__))

Let $\Omega=(\gamma_1,\gamma_2)$ with $\gamma_1=0$, $\gamma_2=1$.
This is the simplest second order boundary value problem (BVP)
for a partial differential equation  (PDE):

```math
-\Delta u =0\\
u(\gamma_1)=g_1\\
u(\gamma_2)=g_2.
```

We replace the Dirichlet boundary condition by a Robin boundary
condition with a penalty parameter $\frac{1}{\varepsilon}$:

```math
\nabla u(\gamma_1) + \frac{1}{\varepsilon}(u(\gamma_1)-g_1)=0  \\
-\nabla u(\gamma_2) + \frac{1}{\varepsilon}(u(\gamma_2)-g_2)
=0  
```

This penalty method for the implementation of Dirichlet
boundary conditions is used throughout VoronoiFVM.

In order to discretize it, we choose collocation points
$\gamma_1=x_1 < x_2 < \dots < x_n=\gamma_2$.

For instance, we can choose 6 collocation points in $(0,1)$:
From these, we create a discretization grid structure
for working with the method.

This implicitly creates a number of control volumes $\omega_k $ around each
discretization point $x_k$: Let
$\sigma_{k,k+1}=\frac{x_k+x_{k+1}}{2}$. Then $\omega_1=(\gamma_1,\sigma_{1,2})$,
$\omega_k= (\sigma_{k-1,k}, \sigma_{k,k+1})$ for $k=2\dots n-1$, $\omega_{n}=(\sigma_{n-1,n},\gamma_2)$. 

```
 x1    x2    x3    x4    x5    x6
 o-----o-----o-----o-----o-----o
 |--|-----|-----|-----|-----|--|
  ω1  ω2     ω3    ω4    ω5  ω6
```

For each $\omega_k$, we integrate the equation

```math
\begin{aligned}
0&=\int_{\omega_k} -\Delta u d\omega=  -\int_{\partial \omega_k} \nabla u ds\\
&= \begin{cases}
u'(\sigma_{1,2}) - u'(0)& k=1\\
u'(\sigma_{k,k+1}) - u'(\sigma_{k-1,k}) & 1<k<n\\
u'(1)- u'(\sigma_{n,n+1})&k=n
\end{cases}\\
&\approx \begin{cases}
\frac{1}{x_2-x_1} g(u_1,u_2) + \frac{1}{\varepsilon}(u_1-0)& k=1\\
\frac{1}{x_k-x_{k-1}}g(u_k,u_{k-1}) -\frac{1}{x_{k+1}-x_{k}}g(u_{k+1},u_{k}) & 1<k<n\\
\frac{1}{\varepsilon}(u_n-1)+ \frac{1}{x_n-x_{n-1}} g(u_{n},u_{n-1})&k=n
\end{cases}
\end{aligned}
```

In the last equation, we wrote $u_k=u(x_k)$ and $g(u_k,u_l)=u_k-u_l$. For
the interior interfaces between control volumes, we replaced $u'$ by a
difference quotient. In the boundary control volumes, we replaced $u'$  by the boundary
conditions.

In the example below, we fix a number of species and  write a Julia function describing $g$,
we create a physics record, and a finite volume system with one unknown species and
a dense matrix to describe it's degrees of freedom
(the matrix used  to calculate the solution is sparse).
We give the species the number 1 and enable it for grid region number one 1.
Then, we set boundary conditions for species 1 at $\gamma_1, \gamma_2$.

We create a zero initial value and a solution vector
and initialize them.

With these data, we solve the system.

We wrap this example and all later ones
into a module structure. This allows to load
all of them at once into the REPL without name
clashes. We shouldn't forget the corresponding end
statement.

=#

module Example101_Laplace1D

using VoronoiFVM, ExtendableGrids

function main()
    ispec = 1    ## Index of species we are working with

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function flux!(f, u, edge, data)
        f[1] = u[1, 1] - u[1, 2]
        return nothing
    end

    function bcond!(f, u, node, data)
        boundary_dirichlet!(f, u, node; region = 1, value = 0)
        boundary_dirichlet!(f, u, node; region = 2, value = 1)
        return nothing
    end

    ## Create a one dimensional discretization grid
    ## Each grid cell belongs to a region marked by a region number
    ## By default, there is only one region numbered with 1
    grid = simplexgrid(0:0.2:1)

    ## Create a finite volume system
    sys = VoronoiFVM.System(grid; flux = flux!, breaction = bcond!, species = ispec)

    ## Solve stationary problem
    solution = solve(sys; inival = 0)

    ## Return test value
    return sum(solution)
end

using Test
function runtests()
    @test main() ≈ 3.0
    return nothing
end

end
