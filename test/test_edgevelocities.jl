module test_edgevelocities
using Test
using LinearAlgebra

using SimplexGridFactory
using Triangulate
using ExtendableGrids

using ExtendableFEM
using ExtendableFEMBase

using VoronoiFVM

function LShape_grid_unsuitable(;
        nref = 1,
        inlet_width = 0.5,
        outlet_width = 0.5,
        inlet_length = 0.5,
        catalyst_length = 0.1,
        outlet_length = 0.5 - catalyst_length,
        maxvol = (2.0)^(-nref - 4),
        refinement_coefficient = 1.0,
        refinement_exponent = 1.0,
        kwargs...
    )

    builder = SimplexGridBuilder(Generator = Triangulate)

    p1 = point!(builder, 0.0, 0.0)
    p2 = point!(builder, inlet_width, 0.0)
    p3 = point!(builder, inlet_width, inlet_length)
    p4c = point!(builder, inlet_width + catalyst_length, inlet_length)
    p4 = point!(builder, inlet_width + catalyst_length + outlet_length, inlet_length)
    p5 = point!(builder, inlet_width + catalyst_length + outlet_length, inlet_length + outlet_width)
    p6 = point!(builder, 0.0, inlet_length + outlet_width)

    facetregion!(builder, 4)
    facet!(builder, p1, p2)

    facetregion!(builder, 2)
    facet!(builder, p4, p5)

    facetregion!(builder, 1)
    facet!(builder, p3, p4c)

    facetregion!(builder, 5)
    facet!(builder, p2, p3)
    facet!(builder, p4c, p4)
    facet!(builder, p5, p6)

    facetregion!(builder, 3)
    facet!(builder, p6, p1)

    points = builder.pointlist.points
    corner_point = points[:, p3]

    function unsuitable(x1, y1, x2, y2, x3, y3, area)
        bary = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3]
        dist_corner = norm(bary - corner_point)
        if area > refinement_coefficient * maxvol * dist_corner^(refinement_exponent)
            return 1
        else
            return 0
        end
    end

    options!(builder, unsuitable = unsuitable)
    grid = simplexgrid(builder)

    return grid
end

function flow_field(x, y)
    return (x^2, -2x * y)
end

function flow_field!(result, qpinfo)
    x = qpinfo.x
    result[1] = x[1]^2
    result[2] = -2 * x[1] * x[2]
    return nothing
end

function test_fem_velocities(; nref = 3, interpolate_eps = 1.0e-10, kwargs...)
    flowgrid = LShape_grid_unsuitable(; nref, kwargs...)

    fespace = FESpace{H1P2{2, 2}}(flowgrid)
    fevector = FEVector(fespace)
    interpolate!(fevector[1], flow_field!)

    evelo_fem = edgevelocities(flowgrid, fevector[1]; interpolate_eps)
    bfvelo_fem = bfacevelocities(flowgrid, fevector[1]; interpolate_eps)

    evelo_analytical = edgevelocities(flowgrid, flow_field)
    bfvelo_analytical = bfacevelocities(flowgrid, flow_field)

    @test norm(evelo_analytical - evelo_fem, Inf) ≤ 1.0e-8
    @test norm(bfvelo_analytical - bfvelo_fem, Inf) ≤ 1.0e-8

    return nothing
end

function runtests()
    test_fem_velocities()
    return nothing
end

end
