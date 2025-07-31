#=

# 311: Heat Equation with boundary diffusion 
([source code](@__SOURCE_URL__))

=#

module Example311_HeatEquation_BoundaryDiffusion
using Printf
using VoronoiFVM
using ExtendableGrids

## Problem data structure to avoid global variables
mutable struct ProblemData
    eps::Float64        # Bulk heat conduction coefficient
    eps_surf::Float64   # Surface diffusion coefficient  
    k::Float64          # Transmission coefficient
    breg::Int           # Boundary region number for surface diffusion
end

"""
  We solve the following system

      ∂_tu - εΔu = 0            in [0,T] × Ω>
           ε∇u⋅ν = k(u-v)       on [0,T] × Γ_1
           ε∇u⋅ν = 0            on [0,T] × (∂Ω ∖ Γ_1)
  ∂_tv -ε_ΓΔ_Γ v = f(x) +k(u-v) on [0,T] × Γ_1
          u(0)   = 0.5          in   {0} × Ω
          v(0)   = 0.5          on   {0} × Γ_1  
"""

function main(n = 1; assembly = :edgewise)
    breg = 5 # boundary region number for surface diffusion

    hmin = 0.05 * 2.0^(-n + 1)
    hmax = 0.2 * 2.0^(-n + 1)
    XLeft = geomspace(0.0, 0.5, hmax, hmin)
    XRight = geomspace(0.5, 1.0, hmin, hmax)
    X = glue(XLeft, XRight)

    Z = geomspace(0.0, 1.0, hmin, 2 * hmax)

    grid = simplexgrid(X, X, Z)

    ## Create problem data structure
    problem_data = ProblemData(1.0e0, 1.0e-2, 1.0, breg)

    physics = VoronoiFVM.Physics(;
        flux = function (f, u, edge, data)
            f[1] = data.eps * (u[1, 1] - u[1, 2])
            return nothing
        end,
        bflux = function (f, u, edge, data)
            if edge.region == data.breg
                f[2] = data.eps_surf * (u[2, 1] - u[2, 2])
            else
                f[2] = 0.0
            end
            return nothing
        end,
        breaction = function (f, u, node, data)
            if node.region == data.breg
                f[1] = data.k * (u[1] - u[2])
                f[2] = data.k * (u[2] - u[1])
            else
                f[1] = 0.0
                f[2] = 0.0
            end
            return nothing
        end,
        bsource = function (f, bnode, data)
            x1 = bnode[1] - 0.5
            x2 = bnode[2] - 0.5
            x3 = bnode[3] - 0.5
            f[2] = 1.0e4 * exp(-20.0 * (x1^2 + x2^2 + x3^2))
            return nothing
        end, bstorage = function (f, u, node, data)
            if node.region == data.breg
                f[2] = u[2]
            end
            return nothing
        end, storage = function (f, u, node, data)
            f[1] = u[1]
            return nothing
        end,
        data = problem_data
    )

    sys = VoronoiFVM.System(grid, physics; unknown_storage = :sparse, assembly)
    enable_species!(sys, 1, [1])
    enable_boundary_species!(sys, 2, [problem_data.breg])

    function tran32!(a, b)
        a[1] = b[2]
        return nothing
    end

    bgrid2 = subgrid(grid, [problem_data.breg]; boundary = true, transform = tran32!)

    U = unknowns(sys)
    U .= 0.5

    control = VoronoiFVM.SolverControl()
    control.verbose = false
    control.reltol_linear = 1.0e-5
    control.keepcurrent_linear = false

    tstep = 0.1
    time = 0.0
    step = 0
    T = 1.0
    while time < T
        time = time + tstep
        U = solve(sys; inival = U, control, tstep)
        tstep *= 1.0
        step += 1
    end

    U_surf = view(U[2, :], bgrid2)
    return sum(U_surf)
end

using Test
function runtests()
    testval = 1509.8109057757858
    testval = 1508.582565216869
    @test isapprox(main(; assembly = :edgewise), testval; rtol = 1.0e-12)
    @test isapprox(main(; assembly = :cellwise), testval; rtol = 1.0e-12)
    return nothing
end

end
