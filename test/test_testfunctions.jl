module test_testfunctions
using VoronoiFVM
using ExtendableGrids
using Test
using LinearAlgebra

function flux!(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

function bcondition!(f, u, bnode, data)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 1, value = 0.0)
    boundary_dirichlet!(f, u, bnode; species = 1, region = 2, value = 1.0)
    return nothing
end

function define_system(n = 5, coordsystem = Cartesian2D, boundary = 2)
    @assert n > 3
    # adaptively refine towards the boundary
    X = geomspace(0.0, 2.0, 2.0^(3 - n), 2.0^(-n))
    Y = linspace(0.0, 2.0, 2^n)

    if boundary == 2
        grid = simplexgrid(X, Y)
    elseif boundary == 3
        grid = simplexgrid(Y, X)
    end

    if coordsystem == Cylindrical2D
        circular_symmetric!(grid)
    end

    system = VoronoiFVM.System(grid; flux = flux!, bcondition = bcondition!, species = [1])
    VoronoiFVM._complete!(system)
    U = unknowns(system; inival = 0.0)

    return system, U
end

function calc_transfer(n = 5, coordsystem = Cartesian2D, boundary = 2; flux_zero = false)
    system, U = define_system(n, coordsystem, boundary)
    grid = system.grid

    testfuncfac = VoronoiFVM.TestFunctionFactory(system)
    tfc_rea = testfunction(testfuncfac, setdiff(1:4, [boundary]), [boundary])
    coords = grid[Coordinates]
    if coordsystem == Cylindrical2D
        tfc_rea .*= coords[1, :]
        cartesian!(grid) # for our proposed method, we need the cartesian node/edgefactors
        system.is_complete = false
        VoronoiFVM._complete!(system) # retrigger computation of node/edgefactors
    end

    if boundary == 2
        if flux_zero
            U[1, :] .= map((x, y) -> (y), grid) # should be 0
        else
            U[1, :] .= map((x, y) -> (x), grid) # should be 1 or 2π
        end
    elseif boundary == 3
        if flux_zero
            # should be 0 for Cylindrical2D, but doesn't work
            U[1, :] .= map((x, y) -> (x), grid)
        else
            U[1, :] .= map((x, y) -> (y), grid) # should be 1 or π
        end
    end

    I = VoronoiFVM.integrate(system, tfc_rea, U)

    if coordsystem == Cylindrical2D
        I .*= 2π
    end

    return I, system, U, tfc_rea
end

function runtests()
    I, _ = calc_transfer(7, Cartesian2D, 2; flux_zero = false)
    @test I[1] ≈ 2.0 atol = 1.0e-1

    I, _ = calc_transfer(7, Cartesian2D, 3; flux_zero = false)
    @test I[1] ≈ 2.0 atol = 1.0e-1

    I, _ = calc_transfer(7, Cylindrical2D, 2; flux_zero = false)
    @test I[1] ≈ 8π atol = 1.0e-1

    I, _ = calc_transfer(7, Cylindrical2D, 3; flux_zero = false)
    @test I[1] ≈ 4π atol = 1.0e-1

    I, _ = calc_transfer(7, Cartesian2D, 2; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-13

    I, _ = calc_transfer(7, Cartesian2D, 3; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-13

    I, _ = calc_transfer(7, Cylindrical2D, 2; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-13

    I, _ = calc_transfer(7, Cylindrical2D, 3; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-1 # why do we lose so much accuracy here compared to the other cases?

    return nothing
end
end
