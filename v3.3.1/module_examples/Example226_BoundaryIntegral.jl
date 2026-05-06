#=

# 226: Terminal flux calculation via test functions, nD, boundary reaction
([source code](@__SOURCE_URL__))

=#

module Example226_BoundaryIntegral

using VoronoiFVM, GridVisualize, ExtendableGrids

function main(;
        n = 10, Plotter = nothing, verbose = false, unknown_storage = :sparse,
        dim = 2, assembly = :edgewise
    )
    n = [101, 21, 5]
    X = collect(range(0.0, 1; length = n[dim]))
    if dim == 1
        grid = simplexgrid(X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [1]
    elseif dim == 2
        grid = simplexgrid(X, X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [4]
    elseif dim == 3
        grid = simplexgrid(X, X, X)
        Γ_where_T_equal_1 = [2]
        Γ_where_T_equal_0 = [4]
    end

    function storage(f, u, node, data)
        f .= u
        return nothing
    end

    function flux(f, u, edge, data)
        f[1] = u[1, 1] - u[1, 2]
        return nothing
    end

    function breaction(f, u, node, data)
        if node.region == Γ_where_T_equal_1[1]
            f[1] = u[1]^2
        end
        return nothing
    end

    physics = VoronoiFVM.Physics(;
        flux = flux,
        storage = storage,
        breaction = breaction
    )

    system = VoronoiFVM.System(grid, physics; assembly = assembly)
    enable_species!(system, 1, [1])
    boundary_dirichlet!(system, 1, Γ_where_T_equal_0[1], 1.0)

    U = solve(system; inival = 0)

    tf = TestFunctionFactory(system)
    T = testfunction(tf, Γ_where_T_equal_0, Γ_where_T_equal_1)

    scalarplot(grid, U[1, :]; Plotter = Plotter, zplane = 0.50001)
    I = integrate(system, T, U)
    B = integrate(system, breaction, U; boundary = true)
    return isapprox(-I[1], B[Γ_where_T_equal_1[1]]; rtol = 1.0e-12)
end

using Test
function runtests()
    @test main(; dim = 1, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 1, unknown_storage = :dense, assembly = :edgewise)
    @test main(; dim = 2, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 2, unknown_storage = :dense, assembly = :edgewise)
    @test main(; dim = 3, unknown_storage = :sparse, assembly = :edgewise)
    @test main(; dim = 3, unknown_storage = :dense, assembly = :edgewise)

    @test main(; dim = 1, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 1, unknown_storage = :dense, assembly = :cellwise)
    @test main(; dim = 2, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 2, unknown_storage = :dense, assembly = :cellwise)
    @test main(; dim = 3, unknown_storage = :sparse, assembly = :cellwise)
    @test main(; dim = 3, unknown_storage = :dense, assembly = :cellwise)
    return nothing
end

end
