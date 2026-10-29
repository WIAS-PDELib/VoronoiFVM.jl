#=

# 424: Initialization of Abstract quantities
 ([source code](@__SOURCE_URL__))

=#
module Example424_AbstractQuantitiesInit
using Printf
using VoronoiFVM
using ExtendableGrids
using GridVisualize
using LinearAlgebra

function main(; N = 5, Plotter = nothing, unknown_storage = :sparse, assembly = :edgewise)
    if 2 * (N ÷ 2) == N
        N = N + 1
    end

    xcoord = range(0, 2; length = N) |> collect
    grid = simplexgrid(xcoord)
    cellmask!(grid, [1], [2], 2)
    system = VoronoiFVM.System(grid; unknown_storage = unknown_storage, assembly = assembly)

    ## First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitly.
    cspec = ContinuousQuantity(system, 1:2)

    ## A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user.
    dspec = DiscontinuousQuantity(system, [1, 2])

    allsubgrids = VoronoiFVM.subgrids(dspec, system)

    function init(u, node)
        ireg = node.region
        return if ireg == 1
            u[dspec] = 1
            u[cspec] = 10
        else
            u[dspec] = 2
            u[cspec] = 20
        end
    end

    function check(u)
        duviews = views(u, dspec, allsubgrids, system)
        cuviews = views(u, cspec, allsubgrids, system)
        result = Bool[]
        psh!(b) = push!(result, b)

        (duviews[1] == fill(1.0, (N - 1) ÷ 2 + 1)) |> psh!
        (duviews[2] == fill(2.0, (N - 1) ÷ 2 + 1)) |> psh!
        (cuviews[2] == fill(20.0, (N - 1) ÷ 2 + 1)) |> psh!
        (cuviews[1][1:(end - 1)] == fill(10.0, (N - 1) ÷ 2)) |> psh!

        return all(result)
    end

    ## "Classical" solution creation
    u = unknowns(system; inifunc = init)

    ## We can use Base.map to create  an initial value
    v = map(init, system)

    ## We also can map  an init function onto the system
    w = unknowns(system)
    map!(init, w, system)

    return check(u) && check(v) && check(w)
end

using Test
function runtests()
    return @test main(; unknown_storage = :sparse, assembly = :edgewise) &&
        main(; unknown_storage = :dense, assembly = :edgewise) &&
        main(; unknown_storage = :sparse, assembly = :cellwise) &&
        main(; unknown_storage = :dense, assembly = :cellwise)
end

end
