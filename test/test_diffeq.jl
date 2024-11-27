module test_diffeq
using VoronoiFVM
using VoronoiFVM: SystemState
using ExtendableGrids
using LinearAlgebra
using Test

function test_matrices(nspec)
    grid = simplexgrid(0:1.0:5)
    function flux(y, u, edge, data)
        for i in 1:length(y)
            y[i] = u[i, 1] - u[i, 2]
        end
        return
    end

    function storage(y, u, node, data)
        for i in 1:length(y)
            y[i] = i * u[i]
        end
        return
    end

    sys = VoronoiFVM.System(grid; flux, storage, species = collect(1:nspec))
    state = SystemState(sys)
    jac_proto = prepare_diffeq!(state, unknowns(sys), 0)
    nd = num_nodes(grid) * nspec
    d = zeros(nd)
    j = 1

    for i in 1:num_nodes(grid)
        fac = 1.0
        if i == 1 || i == num_nodes(grid)
            fac = 0.5
        end
        for id in 1:nspec
            d[j] = fac * id
            j = j + 1
        end
    end
    m = mass_matrix(state)
    u = unknowns(sys)
    J = similar(jac_proto)
    eval_jacobian!(J, u, state, 0.0)
    return @test jac_proto == -J
end

function runtests()
    test_matrices(1)
    test_matrices(2)
    return true
end

end
