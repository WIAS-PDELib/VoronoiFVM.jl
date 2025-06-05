#=

# 205: Convection in axisymmetric stagnation point flow
([source code](@__SOURCE_URL__))

Solve the equation

```math
\begin{aligned}
 -\nabla ( D \nabla u - \vec v u) &= 0\\
             u|_{\Gamma_1} &=1\\
             u|_{\Gamma_0} &=0\\
             (\partial_n u)|_{\Gamma_{out}} & = 0
\end{aligned}
```

in ``\Omega=(0,1)\times (0,1)`` with  ``\Gamma_1 = (0,0.25)\times 1``,
``\Gamma_0=(0.25,1)\times 1`` and  ``\Gamma_{out} = 1\times (0,1)``.
On boundary parts not listed, no-flow boundary conditions are assumed.

The axisymmetric stagnation point flow ``\vec v(r,z)=(vr,-2vz)`` is divergence free.

=#

module Example205_StagnationPoint
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize

function main(; nref = 0, gridname = nothing, Plotter = nothing, D = 0.01, v = 100, cin = 1.0, assembly = :cellwise)
    H = 1.0
    L = 1.0

    Γ_1 = 5
    Γ_0 = 4
    Γ_out = 2

    if !isnothing(gridname)
        grid = simplexgrid(gridname)
    else
        grid = simplexgrid(
            range(0, L; length = 10 * 2^nref),
            range(0, H; length = 10 * 2^nref)
        )
        bfacemask!(grid, [0, H], [0.25L, H], 5)
    end
    circular_symmetric!(grid)

    frz(r, z) = (v * r, -2v * z)

    evelo = edgevelocities(grid, frz)
    bfvelo = bfacevelocities(grid, frz)

    function flux!(f, u, edge, data)
        vd = evelo[edge.index] / D
        bp = fbernoulli(vd)
        bm = fbernoulli(-vd)
        f[1] = D * (bp * u[1] - bm * u[2])
        return nothing
    end

    function outflow!(f, u, node, data)
        if node.region == Γ_out
            f[1] = bfvelo[node.ibnode, node.ibface] * u[1]
        end
        return nothing
    end

    ispec = 1
    physics = VoronoiFVM.Physics(; flux = flux!, breaction = outflow!)
    sys = VoronoiFVM.System(grid, physics; assembly = assembly)
    enable_species!(sys, ispec, [1])
    boundary_dirichlet!(sys, ispec, Γ_1, cin)
    boundary_dirichlet!(sys, ispec, Γ_0, 0)

    tf = TestFunctionFactory(sys)
    tf_in = testfunction(tf, [Γ_out], [Γ_1])
    tf_out = testfunction(tf, [Γ_1], [Γ_out])

    sol = solve(sys)

    I_in = integrate(sys, tf_in, sol)
    I_out = integrate(sys, tf_out, sol)

    scalarplot(sys, sol; Plotter = Plotter)

    ## Test if inflow=outflow
    test1 = isapprox(I_in, -I_out; rtol = 1.0e-5)

    ## Test global maximum principle
    test2 = isapprox(maximum(sol), cin; rtol = 1.0e-10)
    test3 = isapprox(minimum(sol), 0; atol = 1.0e-10)

    ## test zero divergence of fvm velocities
    div = VoronoiFVM.calc_divergences(sys, evelo, bfvelo)
    test4 = all(x -> abs(x) < 1.0e-12, div)

    return test1 && test2 && test3 && test4
end

using Test
function runtests()
    test0 = true
    if VERSION > v"1.6"
        # test on unstructured grid
        gridname = joinpath(pkgdir(VoronoiFVM), "assets", "rz2d.sg")
        test0 = test0 && main(; assembly = :edgewise, gridname) && main(; assembly = :cellwise, gridname)
    end
    @test test0 && main(; assembly = :edgewise) && main(; assembly = :cellwise)
    return nothing
end

end
