#=

# 240: 2D Convection in quadratic stagnation flow velocity field
([source code](@__SOURCE_URL__))

Solve the equation

```math
-\nabla \cdot ( D \nabla u - \mathbf{v} u) = 0
```
in $\Omega=(0,L)\times (0,H)$ with a homogeneous Neumann boundary condition
at $x=0$, an outflow boundary condition at $x=L$, a Dirichlet inflow
condition at $y=H$, and a homogeneous Dirichlet boundary condition
on part of $y=0$.

The analytical expression for the (quadratic variant of the) velocity field 
is $\mathbf{v}(x,y)=(x^2,-2xy)$ in cartesian coordinates and (for the linear variant)
$\mathbf{v}(r,z)=(r,-2z)$ in cylindrical coordinates, i.e.
where the system is solved on $\Omega$ to represent a solution on the solid
of revolution arising from rotating $\Omega$ around $x=0$.

We compute the solution $u$ in both coordinate systems where $\mathbf{v}$ is given
as an analytical expression and as a finite element interpolation onto
the grid of $\Omega$.
=#

module Example240_FiniteElementVelocities
using Printf
using ExtendableFEMBase
using ExtendableFEM
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearAlgebra

function stagnation_flow_cartesian(x, y)
    return (x^2, -2x * y)
end

# In the cylindrical case: since the reconstruction space $\mathtt{HDIVBDM2}$
# is only quadratic, but we have to reconstruct $r \, \mathbf{v}(r,z)$ for a
# (reconstructed) divergence-free solution,
# we can only resolve at most the linear case exactly.
function stagnation_flow_cylindrical(r, z)
    return (r, -2 * z)
end

function inflow_cylindrical(u, qpinfo)
    x = qpinfo.x
    u .= stagnation_flow_cylindrical(x[1], x[2])
    return nothing
end

function inflow_cartesian(u, qpinfo)
    x = qpinfo.x
    u .= stagnation_flow_cartesian(x[1], x[2])
    return nothing
end

function flux!(f, u, edge, data)
    vd = data.evelo[edge.index] / data.D
    bp = fbernoulli(vd)
    bm = fbernoulli(-vd)
    f[1] = data.D * (bp * u[1] - bm * u[2])
    return nothing
end

function bconditions!(f, u, node, data)
    ## catalytic Dirichlet condition at y=0
    if node.region == 5
        boundary_dirichlet!(f, u, node, 1, node.region, 0.0)
    end

    ## outflow condition at x = L
    if node.region == 2
        f[1] = data.bfvelo[node.ibnode, node.ibface] * u[1]
    end

    ## inflow condition at y = H
    if node.region == 3
        boundary_dirichlet!(f, u, node, 1, node.region, data.cin)
    end
    return nothing
end

mutable struct Data
    D::Float64
    cin::Float64
    evelo::Vector{Float64}
    bfvelo::Matrix{Float64}

    Data() = new()
end

#=
Calculate the analytical or FEM solution to the stagnation point flow field $\mathbf{v}$ and
use this as input to solve for the species concentration $u$ of the corresponding
convection-diffusion system.

The passed flags regulate the following behavior:
  - `cylindrical_grid`:     if true, calculates both the velocity field $\mathbf{v}(r,z)$
    and the species concentration $u(r,z)$ in cylindrical coordinates, assuming 
    rotationally symmetric solutions for both.
  - `usefem`:               if true, calculates the velocity field $\mathbf{v}$ using the
    finite element method provided by [`ExtendableFEM`](https://github.com/chmerdon/ExtendableFEM.jl).
  - `reconst`:              if true, interpolates the FEM-calculated 
    velocity field onto a "reconstruction" finite element space that
    provides an exactly divergence-free solution. In a cylindrical grid, 
    this returns not $\mathbf{v}(r,z)$, but $r \, \mathbf{v}(r,z)$ as the velocity.
  - `use_different_grids`:  if true, calculates the FEM solution of the 
    velocity field on a uniform `flowgrid` and the species concentration 
    on an adaptively refined `chemgrid` while still interpolating the 
    calculated velocity correctly onto the `chemgrid`.
=#
function main(;
        cylindrical_grid = false, usefem = true, reconst = cylindrical_grid, use_different_grids = false, nref = 1,
        Plotter = nothing, μ = 1.0e-2, D = 0.01, cin = 1.0, assembly = :edgewise,
        interpolation_eps = 1.0e-9
    )
    H = 1.0
    L = 5.0

    if cylindrical_grid
        coord_system = Cylindrical2D
    else
        coord_system = Cartesian2D
    end

    flowgrid = simplexgrid(
        range(0, L; length = 20 * 2^nref),
        range(0, H; length = 5 * 2^nref)
    )

    if use_different_grids
        h_fine = 1.0e-1
        X_bottom = geomspace(0.0, L / 2, 5.0e-1, h_fine)
        X_cat = range(L / 2, L; step = h_fine)
        chemgrid = simplexgrid(
            [X_bottom; X_cat[2:end]],
            geomspace(0.0, H, 1.0e-3, 1.0e-1)
        )
        bfacemask!(chemgrid, [L / 2, 0.0], [3 * L / 4, 0.0], 5)
    else
        chemgrid = deepcopy(flowgrid)
        bfacemask!(chemgrid, [L / 2, 0.0], [3 * L / 4, 0.0], 5)
    end

    if usefem
        velocity = compute_velocity(
            flowgrid, cylindrical_grid, reconst, μ; interpolation_eps
        )
        DivIntegrator = L2NormIntegrator([div(1)]; quadorder = 2 * 2, resultdim = 1)
        div_v = sqrt(sum(evaluate(DivIntegrator, [velocity])))
        @info "||div(R(v))||_2 = $(div_v)"
    else
        if cylindrical_grid
            velocity = stagnation_flow_cylindrical
        else
            velocity = stagnation_flow_cartesian
        end
    end

    if cylindrical_grid
        analytical_velocity = stagnation_flow_cylindrical
    else
        analytical_velocity = stagnation_flow_cartesian
    end

    ## only the chemgrid needs its CoordinateSystem adjusted
    ## since the velocity calculation works by adjusting
    ## the kernels for the Stokes operator directly while
    ## the finite volume machinery relies upon the CoordinateSystem
    ## for selecting the correct geometrical factors for the
    ## Voronoi cell contributions
    if cylindrical_grid
        chemgrid[CoordinateSystem] = Cylindrical2D
    end

    data = Data()
    data.D = D
    data.cin = cin

    evelo = edgevelocities(chemgrid, velocity; reconst)
    bfvelo = bfacevelocities(chemgrid, velocity; reconst)

    data.evelo = evelo
    data.bfvelo = bfvelo

    physics = VoronoiFVM.Physics(; flux = flux!, breaction = bconditions!, data)
    sys = VoronoiFVM.System(chemgrid, physics; assembly)
    enable_species!(sys, 1, [1])

    sol = solve(sys; inival = 0.0)

    fvm_divs = VoronoiFVM.calc_divergences(sys, evelo, bfvelo)
    @info "||div(v)||_∞ = $(norm(fvm_divs, Inf))"

    vis = GridVisualizer(; Plotter = Plotter)

    scalarplot!(
        vis[1, 1], chemgrid, sol[1, :]; flimits = (0, cin + 1.0e-5),
        show = true
    )

    minmax = extrema(sol)
    @info "Minimal/maximal values of concentration: $(minmax)"

    return sol, evelo, bfvelo, minmax
end

using Test
function runtests()
    cin = 1.0
    for cylindrical_grid in [false, true]
        sol_analytical, evelo_analytical, bfvelo_analytical, minmax_analytical = main(;
            cylindrical_grid, cin, usefem = false
        )
        sol_fem, evelo_fem, bfvelo_fem, minmax_fem = main(;
            cylindrical_grid, cin, usefem = true
        )
        @test norm(evelo_analytical .- evelo_fem, Inf) ≤ 1.0e-11
        @test norm(bfvelo_analytical .- bfvelo_fem, Inf) ≤ 1.0e-9
        @test norm(sol_analytical .- sol_fem, Inf) ≤ 1.0e-10
        @test norm(minmax_analytical .- [0.0, cin], Inf) ≤ 1.0e-15
        @test norm(minmax_fem .- [0.0, cin], Inf) ≤ 1.0e-11
    end
    return nothing
end

function compute_velocity(
        flowgrid, cylindrical_grid, reconst, μ = 1.0e-2; interpolation_eps = 1.0e-10
    )
    ## define finite element spaces
    FE_v, FE_p = H1P2B{2, 2}, L2P1{1}
    reconst_FEType = HDIVBDM2{2}
    FES = [FESpace{FE_v}(flowgrid), FESpace{FE_p}(flowgrid; broken = true)]

    ## describe problem
    Problem = ProblemDescription("incompressible Stokes problem")
    v = Unknown("v"; name = "velocity")
    p = Unknown("p"; name = "pressure")
    assign_unknown!(Problem, v)
    assign_unknown!(Problem, p)

    ## assign stokes operator
    assign_operator!(
        Problem,
        BilinearOperator(
            kernel_stokes!, cylindrical_grid ? [id(v), grad(v), id(p)] : [grad(v), id(p)];
            bonus_quadorder = 2, store = false,
            params = [μ, cylindrical_grid]
        )
    )

    ## assign Dirichlet boundary conditions on all boundary regions to
    ## enforce match with analytical solution
    if cylindrical_grid
        assign_operator!(
            Problem, InterpolateBoundaryData(v, inflow_cylindrical; regions = [1, 2, 3, 4])
        )
    else
        assign_operator!(
            Problem, InterpolateBoundaryData(v, inflow_cartesian; regions = [1, 2, 3, 4])
        )
    end

    velocity_solution = solve(Problem, FES)

    ## ensure divergence free solution by projecting onto reconstruction spaces
    FES_reconst = FESpace{reconst_FEType}(flowgrid)
    R = FEVector(FES_reconst)
    if reconst
        if cylindrical_grid
            lazy_interpolate!(
                R[1], velocity_solution, [id(v)]; postprocess = multiply_r,
                bonus_quadorder = 2, eps = interpolation_eps
            )
        else
            lazy_interpolate!(
                R[1], velocity_solution, [id(v)];
                bonus_quadorder = 2, eps = interpolation_eps
            )
        end
    else
        return velocity_solution[1]
    end

    return R[1]
end

function kernel_stokes!(result, u_ops, qpinfo)
    μ = qpinfo.params[1]
    cylindrical_grid = qpinfo.params[2]
    if cylindrical_grid > 0
        r = qpinfo.x[1]
        u, ∇u, p = view(u_ops, 1:2), view(u_ops, 3:6), view(u_ops, 7)
        result[1] = μ / r * u[1] - p[1]
        result[2] = 0
        result[3] = μ * r * ∇u[1] - r * p[1]
        result[4] = μ * r * ∇u[2]
        result[5] = μ * r * ∇u[3]
        result[6] = μ * r * ∇u[4] - r * p[1]
        result[7] = -(r * (∇u[1] + ∇u[4]) + u[1])
    else
        ∇u, p = view(u_ops, 1:4), view(u_ops, 5)
        result[1] = μ * ∇u[1] - p[1]
        result[2] = μ * ∇u[2]
        result[3] = μ * ∇u[3]
        result[4] = μ * ∇u[4] - p[1]
        result[5] = -(∇u[1] + ∇u[4])
    end
    return nothing
end

function multiply_r(result, input, qpinfo)
    x = qpinfo.x
    result .= input * x[1]
    return nothing
end

end
