#=

# 203: Various coordinate systems 
([source code](@__SOURCE_URL__))
=#

module Example203_CoordinateSystems

using VoronoiFVM
using LinearAlgebra
using ExtendableGrids
using GridVisualize

function plot(grid, numerical, exact, Plotter)
    vis = GridVisualizer(; Plotter = Plotter, layout = (2, 1))
    scalarplot!(vis[1, 1], grid, numerical[1, :]; title = "numerical")
    return scalarplot!(vis[2, 1], grid, exact; title = "exact", show = true)
end

function flux(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

"""
    symlapdisk(r,r2)

Exact solution of homogeneous Dirichlet problem `-Δu=1` on disk of radius r2. 
"""
symlapdisk(r, r2) = 0.25 * (r2^2 - r^2)

"""
    maindisk(;nref=0, r2=5.0, Plotter=nothing)

Solve homogeneuous Dirichlet problem  `-Δu=1`  
on disk of radius r2, exact solution is `(r_2^2-r^2)/4`. 

In this case, the discretization appears to be exact.
"""
function maindisk(; nref = 0, r2 = 5.0, Plotter = nothing, assembly = :edgewise)
    h = 0.1 * 2.0^(-nref)
    R = collect(0:h:r2)
    grid = simplexgrid(R)
    circular_symmetric!(grid)
    source(f, node, data) = f[1] = 1.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapdisk.(coordinates(grid)[1, :], r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) < 1.0e-14
end

"""
    maincylinder(;nref=0, r2=5.0, z1=0, z2=1, Plotter=nothing)

Solve homogeneuous Dirichlet problem  `-Δu=1`  
on disk of radius r2, exact solution is `(r_2^2-r^2)/4`. 

In this case, the discretization appears to be exact.
"""
function maincylinder(;
        nref = 0,
        r2 = 5.0,
        z1 = 0.0,
        z2 = 1.0,
        Plotter = nothing,
        assembly = :edgewise,
    )
    h = 0.1 * 2.0^(-nref)
    R = collect(0:h:r2)
    Z = collect(z1:h:z2)
    grid = simplexgrid(R, Z)
    circular_symmetric!(grid)
    source(f, node, data) = f[1] = 1.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapdisk.(coordinates(grid)[1, :], r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) < 1.0e-14
end

"""
    maincylinder_unstruct(;Plotter=nothing)

Solve homogeneuous Dirichlet problem  `-Δu=1`  
on disk of radius r2, exact solution is `(r_2^2-r^2)/4`. 

In this case, the discretization appears to be exact.
"""
function maincylinder_unstruct(;
        Plotter = nothing,
        assembly = :edgewise
    )
    if VERSION < v"1.7"
        # no pkdir
        return true
    end
    nref = 0
    r2 = 5.0
    z1 = 0.0
    z2 = 1.0
    h = 0.1 * 2.0^(-nref)
    grid = simplexgrid(joinpath(pkgdir(VoronoiFVM), "assets", "cyl_unstruct.sg"))
    circular_symmetric!(grid)
    source(f, node, data) = f[1] = 1.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapdisk.(coordinates(grid)[1, :], r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) < 0.0012
end

"""
    symlapring(r,r1,r2)

Exact solution of Dirichlet problem `-Δu=0` on ring between radii r1 and r2, 
with boundary value 1 at r1 and 0 at r2. 
"""
symlapring(r, r1, r2) = (log(r) - log(r2)) / (log(r1) - log(r2))

"""
    mainring(;nref=0, r2=5.0, Plotter=nothing)

of Dirichlet problem `-Δu=0` on ring between radii r1 and r2, 
with boundary value 1 at r1 and 0 at r2. Test of quadratic convergence.
"""
function mainring(; nref = 0, r1 = 1.0, r2 = 5.0, Plotter = nothing, assembly = :edgewise)
    h = 0.1 * 2.0^(-nref)
    R = collect(r1:h:r2)
    grid = simplexgrid(R)
    circular_symmetric!(grid)
    source(f, node, data) = f[1] = 0.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 1, value = 1.0)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapring.(coordinates(grid)[1, :], r1, r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) / h^2 < 0.01
end

"""
    maincylindershell(;nref=0, r2=5.0, z1=0.0, z2=1.0, Plotter=nothing)

of Dirichlet problem `-Δu=0` on cylindershell between radii r1 and r2, 
with boundary value 1 at r1 and 0 at r2. Test of quadratic convergence.
"""
function maincylindershell(;
        nref = 0,
        r1 = 1.0,
        r2 = 5.0,
        z1 = 0.0,
        z2 = 1.0,
        Plotter = nothing,
        assembly = :edgewise,
    )
    h = 0.1 * 2.0^(-nref)
    R = collect(r1:h:r2)
    Z = collect(z1:h:z2)
    grid = simplexgrid(R, Z)
    circular_symmetric!(grid)
    source(f, node, data) = f[1] = 0.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 4, value = 1.0)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapring.(coordinates(grid)[1, :], r1, r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) / h^2 < 0.01
end

"""
    symlapsphere(r,r2)

Exact solution of homogeneous Dirichlet problem `-Δu=1` on sphere of radius r2. 
"""
symlapsphere(r, r2) = (r2^2 - r^2) / 6.0

"""
    mainsphere(;nref=0, r2=5.0, Plotter=nothing)

Solve homogeneuous Dirichlet problem  `-Δu=1`  
on sphere of radius r2, exact solution is `(r_2^2-r^2)/4`. 

In this case, the discretization appears to be exact.
"""
function mainsphere(; nref = 0, r2 = 5.0, Plotter = nothing, assembly = :edgewise)
    h = 0.1 * 2.0^(-nref)
    R = collect(0:h:r2)
    grid = simplexgrid(R)
    spherical_symmetric!(grid)
    source(f, node, data) = f[1] = 1.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapsphere.(coordinates(grid)[1, :], r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) < 1.0e-14
end

"""
    symlapsphereshell(r,r1,r2)

Exact solution of Dirichlet problem `-Δu=0` on sphereshell between radii r1 and r2, 
with boundary value 1 at r1 and 0 at r2. 
"""
symlapsphereshell(r, r1, r2) = (r2 * r1 / r - r1) / (r2 - r1)

"""
    mainsphereshell(;nref=0, r2=5.0, Plotter=nothing)

of Dirichlet problem `-Δu=0` on sphereshell between radii r1 and r2, 
with boundary value 1 at r1 and 0 at r2. Test of quadratic convergence.
"""
function mainsphereshell(;
        nref = 0,
        r1 = 1.0,
        r2 = 5.0,
        Plotter = nothing,
        assembly = :edgewise,
    )
    h = 0.1 * 2.0^(-nref)
    R = collect(r1:h:r2)
    grid = simplexgrid(R)
    spherical_symmetric!(grid)
    source(f, node, data) = f[1] = 0.0
    sys = VoronoiFVM.System(grid; source, flux, species = [1], assembly = assembly)
    boundary_dirichlet!(sys; species = 1, region = 1, value = 1.0)
    boundary_dirichlet!(sys; species = 1, region = 2, value = 0.0)
    sol = solve(sys)
    exact = symlapsphereshell.(coordinates(grid)[1, :], r1, r2)
    plot(grid, sol, exact, Plotter)
    return norm(sol[1, :] - exact, Inf) / h^2 < 0.04
end

#
# Called by unit test

using Test #
function runtests()
    @test maindisk(; assembly = :edgewise) &&
        mainring(; assembly = :edgewise) &&
        maincylinder(; assembly = :edgewise) &&
        maincylinder_unstruct(; assembly = :edgewise) &&
        maincylindershell(; assembly = :edgewise) &&
        mainsphere(; assembly = :edgewise) &&
        mainsphereshell(; assembly = :edgewise) &&
        maindisk(; assembly = :cellwise) &&
        mainring(; assembly = :cellwise) &&
        maincylinder(; assembly = :cellwise) &&
        maincylinder_unstruct(; assembly = :cellwise) &&
        maincylindershell(; assembly = :cellwise) &&
        mainsphere(; assembly = :cellwise) &&
        mainsphereshell(; assembly = :cellwise)
    return nothing
end

end
