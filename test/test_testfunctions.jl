module test_testfunctions
using VoronoiFVM
using ExtendableGrids
using Test
using LinearAlgebra

function flux!(f, u, edge, data)
    f[1] = u[1, 1] - u[1, 2]
    return nothing
end

function _source_closure(func)
    return function source!(f, node, data)
        return f[1] = func(node[:]...)
    end
end

function source_closure(coordsystem = Cartesian2D, boundary = 2, flux_zero = false)
    if coordsystem == Cartesian2D
        return _source_closure((x, y) -> (0))
    elseif coordsystem == Cylindrical2D
        if boundary == 2
            if flux_zero
                #return _source_closure((r, z) -> (-2 * z + 2 * z / (r + 1.0e-9)))
                return _source_closure((r, z) -> (0.0))
            else
                return _source_closure((r, z) -> 0.0)
            end
        elseif boundary == 3
            if flux_zero
                #return _source_closure((r, z) -> (-r - z^2 / (2 * (r + 1.0e-8)) + 2 * z / (r + 1.0e-8)))
                return _source_closure((r, z) -> (-1 / (r + 1.0e-8)))
            else
                return _source_closure((r, z) -> 0.0)
            end
        end
    end
end

function solution_closure(coordsystem = Cartesian2D, boundary = 2, flux_zero = false)
    return if coordsystem == Cartesian2D
        if boundary == 2
            if flux_zero
                return (x, y) -> y
            else
                return (x, y) -> x
            end
        elseif boundary == 3
            if flux_zero
                return (x, y) -> x
            else
                return (x, y) -> y
            end
        end
    elseif coordsystem == Cylindrical2D
        if boundary == 2
            if flux_zero
                #return (r, z) -> ((z / 2) * r^2 - 2 * r * z)
                return (r, z) -> (z)
            else
                return (r, z) -> (0.5 * r^2 - z^2)
            end
        elseif boundary == 3
            if flux_zero
                #return (r, z) -> ((r / 2) * z^2 - 2 * r * z)
                return (r, z) -> (r)
            else
                return (r, z) -> (0.5 * r^2 - z^2)
            end
        end
    end
end

function define_system(n = 5, coordsystem = Cartesian2D, boundary = 2; sol, flux_zero = false)
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

    src! = source_closure(coordsystem, boundary, flux_zero)

    system = VoronoiFVM.System(grid; flux = flux!, source = src!, species = [1])
    VoronoiFVM._complete!(system)
    U = unknowns(system; inival = 0.0)

    return system, U
end

function calc_transfer(n = 5, coordsystem = Cartesian2D, boundary = 2; flux_zero = false)
    # calculate ansatz function for given case
    sol = solution_closure(coordsystem, boundary, flux_zero)

    # construct system and interpolate ansatz function onto grid nodes
    system, U = define_system(n, coordsystem, boundary; sol, flux_zero)
    grid = system.grid
    U[1, :] .= map(sol, grid)

    # generate test function
    testfuncfac = VoronoiFVM.TestFunctionFactory(system)
    tfc_rea = testfunction(testfuncfac, setdiff(1:4, [boundary]), [boundary])

    I = VoronoiFVM.integrate(system, tfc_rea, U)

    return I, system, U, tfc_rea
end

function runtests()
    I, _ = calc_transfer(7, Cartesian2D, 2; flux_zero = false)
    @test I[1] ≈ 2.0 atol = 1.0e-1

    I, _ = calc_transfer(7, Cartesian2D, 3; flux_zero = false)
    @test I[1] ≈ 2.0 atol = 1.0e-1

    I, _ = calc_transfer(7, Cylindrical2D, 2; flux_zero = false)
    @test I[1] ≈ 16π atol = 5.0e-1

    I, _ = calc_transfer(7, Cylindrical2D, 3; flux_zero = false)
    @test I[1] ≈ -16π atol = 5.0e-1

    I, _ = calc_transfer(7, Cartesian2D, 2; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-13

    I, _ = calc_transfer(7, Cartesian2D, 3; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-13

    I, _ = calc_transfer(7, Cylindrical2D, 2; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-1 # why does the accuracy decrease so much?

    I, _ = calc_transfer(7, Cylindrical2D, 3; flux_zero = true)
    @test I[1] ≈ 0.0 atol = 1.0e-1 # why does this test case fail?

    return nothing
end
end
