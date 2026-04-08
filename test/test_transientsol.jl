module test_transientsol
using VoronoiFVM
using SparseArrays
using Test
function make_transientsol(; n = 10, M = 5, N = 100, in_memory = true)
    makevec(k) = [k + i * j for i in 1:M, j in 1:N]
    sol = TransientSolution(0, makevec(0); in_memory = in_memory)
    for k in 1:n
        append!(sol, k, makevec(k))
    end
    return sol
end

function make_transientsol2(;
        n = 10, M = 5, N = 100,
        in_memory = true,
        sparsestorage::Bool = false,
        keep_open::Bool = false
    )
    function makevec(k)
        v = Float64[k + i * j for i in 1:M, j in 1:N]
        if sparsestorage
            v = sparse(v)
        end
        return v
    end

    sol = TransientSolution(0.0, makevec(0); in_memory, keep_open)
    for k in 1:n
        append!(sol, k, makevec(k))
    end
    return sol
end


function runtests()
    msol = make_transientsol(; in_memory = true)
    dsol = make_transientsol(; in_memory = false)
    @test msol == dsol
    @test length(msol[2, :, 1]) == 100
    @test length(dsol[2, :, 2]) == 100

    for keep_open in [true, false]
        for sparsestorage in [false, true]
            msol = make_transientsol2(; in_memory = true, sparsestorage, keep_open)
            dsol = make_transientsol2(; in_memory = false, sparsestorage, keep_open)
            @test msol == dsol
            for sol in [msol, dsol]
                @test size(sol.u[1]) == size(sol.u[10])
                @test size(sol.u[1]) == (5, 100)
                @test size(sol[:, :, 5]) == (5, 100)
                @test size(sol[3, :, :]) == (100, 11)
                @test size(sol[:, 7, :]) == (5, 11)
                @test sol[:, 3, 4] == 6:3:18
                @test sol[3, :, 4] == 6:3:303
                @test sol[2, 5, :] == 10:20
                @test sol(3.5) ≈ sol[:, :, 3] .+ 1.5
                @test length(sol[2, :, 1]) == 100
            end
        end
    end
    return true
end

end
