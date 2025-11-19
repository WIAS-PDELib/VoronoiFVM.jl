#=

# 410: Many Species
 ([source code](@__SOURCE_URL__))

Test stationary diffusion for 50 species.

=#
module Example410_ManySpecies
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids
using GridVisualize
using LinearAlgebra

function main(; n = 11, nspec = 50, Plotter = nothing, unknown_storage = :dense, assembly = :edgewise)
    grid = simplexgrid(range(0, 1; length = n))

    function flux(f, u, edge, data)
        for ispec in 1:nspec
            f[ispec] = u[ispec, 1] - u[ispec, 2]
        end
        return nothing
    end
    physics = VoronoiFVM.Physics(; flux = flux)
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)
    for ispec in 1:nspec
        enable_species!(sys, ispec, [1])
        boundary_dirichlet!(sys, ispec, 1, 0)
        boundary_dirichlet!(sys, ispec, 2, 1)
    end
    sol = solve(sys)
    return norm(sol)
end

using Test
function runtests()
    testval = 13.874436925511608
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
