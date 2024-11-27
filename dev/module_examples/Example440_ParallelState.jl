#=

# 440: Parallel solves
 ([source code](@__SOURCE_URL__))

Demonstrate how to solve one system with different data in parallel using the SystemState (new in v2.0).
=#

module Example440_ParallelState


using VoronoiFVM, ExtendableGrids
using GridVisualize
using ChunkSplitters

function flux(y, u, node, data)
    y[1] = u[1, 1]^2 - u[1, 2]^2
    return nothing
end

function bcondition(y, u, node, data)
    boundary_dirichlet!(y, u, node, species = 1, region = 1, value = 0.1)
    boundary_neumann!(y, u, node, species = 1, region = 2, value = data.influx)
    return nothing
end

function main(; nref = 5, Plotter = nothing)
    grid = simplexgrid(range(0, 1, length = 10 * 2^nref + 1))
    sys = VoronoiFVM.System(grid; flux, bcondition, species = [1], data = (influx = 0.0,))

    ## Initial state. First solution creates the matrix
    state0 = VoronoiFVM.SystemState(sys)
    sol = solve!(state0; inival = 0.1)

    ## Prepare parameter and result data
    influxes = range(0.0, 10.0, length = 100)
    masses = similar(influxes)

    ## Split the index range in as many chunks as threads
    Threads.@threads for indexes in chunks(1:length(influxes); n = Threads.nthreads())
        ## Create a new state sharing the system - one for each chunk
        state = similar(state0)
        ## Solve for all data values in chunk
        for iflux in indexes
            data = (influx = influxes[iflux],)
            sol = solve!(state; data, inival = 0.1, verbose = "")
            masses[iflux] = integrate(sys, sol)[1, 1]
        end
    end
    scalarplot(influxes, masses; Plotter, xlabel = "influx", ylabel = "mass")
    return sum(masses)
end

using Test
function runtests()
    testval = 140.79872772042577
    @test main() ≈ testval
    return nothing
end

end
