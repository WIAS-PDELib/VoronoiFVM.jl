module test_state
using VoronoiFVM
using VoronoiFVM: SystemState
using ExtendableGrids
using LinearAlgebra
using Test

flux(y, u, edge, data) = y[1] = u[1, 1] - u[1, 2]

function bcondition(y, u, bnode, data)
    boundary_robin!(y, u, bnode, region = 1, value = 0.0, factor = 0.1)
    return boundary_robin!(y, u, bnode, region = 2, value = 1.0, factor = 0.1)
end

function main(; unknown_storage = :dense)
    g = simplexgrid(0:0.1:1)
    sys = VoronoiFVM.System(g; flux, bcondition, species = [1], unknown_storage)
    sol1 = solve(sys)

    # Solution and solution with state shall be the same
    state = VoronoiFVM.SystemState(sys)
    sol2 = solve!(state)
    @test sol1 ≈ sol2

    control = SolverControl()
    fixed_timesteps!(control, 0.025)


    # Allow initial values as result of previous time evolution
    tsol1 = solve(sys; inival = 0.0, times = (0, 0.1), control)
    tsol2 = solve(sys; inival = tsol1.u[end], times = (0.1, 0.2), control)
    tsol3 = solve(sys; inival = 0.0, times = [0.0, 0.1, 0.2], control)
    @test tsol3.u[end] ≈ tsol2.u[end]

    # Solution and solution with state shall be the same
    xsol1 = solve!(state; inival = 0.0, times = (0, 0.1), control)
    xsol2 = solve!(state; inival = tsol1.u[end], times = (0.1, 0.2), control)
    xsol3 = solve!(state; inival = 0.0, times = [0.0, 0.1, 0.2], control)
    @test xsol3.u[end] ≈ xsol2.u[end]
    return @test xsol3.u[end] ≈ tsol3.u[end]


end

function runtests()
    main(; unknown_storage = :dense)
    return main(; unknown_storage = :sparse)

end
end
