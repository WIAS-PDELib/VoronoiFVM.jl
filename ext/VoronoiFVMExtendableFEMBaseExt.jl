module VoronoiFVMExtendableFEMBaseExt

using StaticArrays: @MVector
using VoronoiFVM
using ExtendableFEMBase: CellFinder, evaluate_bary!,
    FEVectorBlock, initialize!, Identity,
    PointEvaluator, SegmentIntegrator

using ExtendableGrids: ExtendableGrid, Cartesian2D, CoordinateSystem, Cylindrical2D,
    Edge1D, eval_trafo!, gFindLocal!,
    mapderiv!, postprocess_xreftest!, update_trafo!

using LinearAlgebra: dot

using Base: fill!

using DocStringExtensions: DocStringExtensions, SIGNATURES

id(u) = (u, Identity)


dist(p1, p2) = sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)

function iscloser(pint, p1, p2, eps)
    return dist(pint, p2) < dist(p2, p1) - eps
end

#This is the FEVectorBlock with lots of info from ExtendableFEMBase
#to execute integration along given segments
struct AugmentedFEVectorBlock{TVB, TSE, TPE, TCF, TFG}
    vblock::TVB
    seg_integrator::TSE
    point_evaluator::TPE
    cellfinder::TCF
    flowgrid::TFG

    bp1::Vector{Float64}
    bp2::Vector{Float64}
    result::Vector{Float64}
    summand::Vector{Float64}
    pint::Vector{Float64}
    bpint::Vector{Float64}
    bary::Vector{Float64}


end

function AugmentedFEVectorBlock(vel, seg_integrator, point_evaluator, cf, flowgrid)


    bp1 = zeros(Float64, 3)
    bp2 = zeros(Float64, 3)
    result = zeros(Float64, 2)
    summand = zeros(Float64, 2)
    pint = zeros(Float64, 2)
    bpint = zeros(Float64, 3)
    bary = [1 / 3, 1 / 3, 1 / 3]


    return AugmentedFEVectorBlock(
        vel, seg_integrator, point_evaluator, cf, flowgrid,
        bp1,
        bp2,
        result,
        summand,
        pint,
        bpint,
        bary
    )


end

function multiply_r(result, input, qpinfo)
    x = qpinfo.x
    result .= input * x[1]
    return nothing
end

function prepare_segment_integration(vel; axisymmetric = false, reconst = false, kwargs...)
    flowgrid = vel.FES.xgrid

    # reference mappings not implemented for other coord types
    cartesian!(flowgrid)

    if axisymmetric
        if reconst
            seg_integrator = SegmentIntegrator(Edge1D, [id(1)])
        else
            seg_integrator = SegmentIntegrator(Edge1D, multiply_r, [id(1)]; bonus_quadorder = 1)
        end
    else
        seg_integrator = SegmentIntegrator(Edge1D, [id(1)])
    end

    initialize!(seg_integrator, [vel])
    point_evaluator = PointEvaluator([id(1)], [vel])

    cellfinder = CellFinder(flowgrid)

    if axisymmetric
        circular_symmetric!(flowgrid)
    end

    return seg_integrator, point_evaluator, cellfinder, flowgrid
end

"""
$(SIGNATURES)

Compute [`VoronoiFVM.edgevelocities`](@ref) for a finite element flow field computed
by [`ExtendableFEM`](https://github.com/chmerdon/ExtendableFEM.jl).
"""
function VoronoiFVM.edgevelocities(grid::ExtendableGrid, vel::FEVectorBlock; kwargs...)
    # construct an augmented type to gather precomputed information
    # in order to pass it to the repeated integrate call in VoronoiFVM.edgevelocities
    axisymmetric = grid[CoordinateSystem] <: Cylindrical2D ? true : false
    seg_integrator, point_evaluator, cf, flowgrid = prepare_segment_integration(vel; axisymmetric, kwargs...)
    aug_fevec_block = AugmentedFEVectorBlock(vel, seg_integrator, point_evaluator, cf, flowgrid)
    velovec = VoronoiFVM.edgevelocities(grid, aug_fevec_block; axisymmetric, kwargs...)

    if axisymmetric
        circular_symmetric!(flowgrid)
    end

    return velovec
end

"""
$(SIGNATURES)

Compute [`VoronoiFVM.bfacevelocities`](@ref) for a finite element flow field computed
by [`ExtendableFEM`](https://github.com/chmerdon/ExtendableFEM.jl).
"""
function VoronoiFVM.bfacevelocities(grid::ExtendableGrid, vel::FEVectorBlock; kwargs...)
    axisymmetric = grid[CoordinateSystem] <: Cylindrical2D ? true : false
    seg_integrator, point_evaluator, cf, flowgrid = prepare_segment_integration(vel; axisymmetric, kwargs...)
    aug_fevec_block = AugmentedFEVectorBlock(vel, seg_integrator, point_evaluator, cf, flowgrid)

    velovec = VoronoiFVM.bfacevelocities(grid, aug_fevec_block; axisymmetric, kwargs...)

    if axisymmetric
        circular_symmetric!(flowgrid)
    end

    return velovec
end

# We need two explicitly type-annotated methods for a working method specialization.
# This one...
function VoronoiFVM.integrate(::Type{<:Cartesian2D}, p1, p2, hnormal, aug_vec_block::AugmentedFEVectorBlock; kwargs...)
    return _integrate_along_segments(p1, p2, hnormal, aug_vec_block; kwargs...)
end

# ... and that one.
function VoronoiFVM.integrate(::Type{<:Cylindrical2D}, p1, p2, hnormal, aug_vec_block::AugmentedFEVectorBlock; kwargs...)
    return _integrate_along_segments(p1, p2, hnormal, aug_vec_block; kwargs...)
end

# compute the path integral for the velocity in aug_vec_block between p1 and p2 by
# incrementally walking through each cell in the grid between p1 and p2
# and summing up each cell's contribution
function _integrate_along_segments(p1, p2, hnormal, aug_vec_block::AugmentedFEVectorBlock{TVB, TSE, TPE, TCF, TFG}; interpolate_eps = 1.0e-12, axisymmetric = false, kwargs...) where {TVB, TSE, TPE, TCF, TFG}
    edge_length = dist(p1, p2)
    avg_r = (p1[1] + p2[1]) / 2
    (; bp1, result, summand, bp2, pint, bpint, bary) = aug_vec_block


    if axisymmetric && avg_r < eps()
        return 0
    end

    CF = aug_vec_block.cellfinder
    icell::Int = gFindLocal!(bp1, CF, p1; eps = interpolate_eps)
    if edge_length ≤ interpolate_eps
        point_evaluator = aug_vec_block.point_evaluator
        evaluate_bary!(p2, point_evaluator, bp1, icell)
        if axisymmetric
            return dot(p2, hnormal) / (avg_r)
        else
            return dot(p2, hnormal)
        end
    end

    SI = aug_vec_block.seg_integrator


    xCellFaces::Matrix{Int} = CF.xCellFaces
    xFaceCells::Matrix{Int} = CF.xFaceCells
    facetogo = CF.facetogo
    cx = CF.cx
    icell_new = icell
    prevcell = icell
    L2G = CF.L2G4EG[1]
    L2Gb::Vector{Float64} = L2G.b
    invA = CF.invA

    t = 0.0

    p1_temp2 = @MVector zeros(2)
    p1_temp3 = @MVector zeros(3)

    function calc_barycentric_coords!(bp, p)
        for j in 1:2
            cx[j] = p[j] - L2Gb[j]
        end

        fill!(bp, 0)
        for j in 1:2, k in 1:2
            bp[k] += invA[j, k] * cx[j]
        end
        return postprocess_xreftest!(bp, CF.xCellGeometries[icell])
    end

    while (true)

        # TODO implement proper emergency guard to avoid indefinite loops

        # first compute the barycentric coordinates of
        # p1,p2

        # update local 2 global map
        L2G = CF.L2G4EG[1]
        update_trafo!(L2G, icell)
        L2Gb = L2G.b

        mapderiv!(invA, L2G, p1)

        calc_barycentric_coords!(bp1, p1)
        calc_barycentric_coords!(bp2, p2)

        # if p1 is a node of a triangle, start with
        # a cell containing p1 in the direction of (p2-p1)
        if count(<=(interpolate_eps), bp1) == 2 # 10^(-13)
            @. p1_temp2 = p1 + 10 * interpolate_eps * (p2 - p1)
            icell_new = gFindLocal!(bp1, CF, p1_temp2; eps = 10 * interpolate_eps, icellstart = icell) #!!! allocates
            if icell_new == 0
                # TODO: test the following
                # icell_new = gFindBruteForce!(bp1, CF, p1_temp[1:2])
                @warn "icell_new=0!"
            end
            if icell_new != icell
                icell = icell_new
                continue
            end
        end

        # push p1 a little towards the triangle circumcenter
        # to avoid it being situated on an edge or node of the
        # flowgrid
        # (to avoid being stuck on an edge if (p2-p1) is
        # on the edge
        @. bp1 += interpolate_eps * (bary - bp1)
        eval_trafo!(p1, L2G, bp1)

        (λp2min, imin) = findmin(bp2)

        # if λp2min≥0, then p2 is inside icell and we can simply add the
        # integral across the line segment between p1 and p2 to the result

        # if not, then p2 is outside of icell and we try to determine
        # pint which is the point where icell intersects the line segment
        # [p1,p2] - since pint should be on the boundary of icell,
        # at least one barycentric coordinate (stored in bpint) should
        # be zero yielding an expression for the line segment parameter t.
        # this is not necessarily the previous imin and we have to check
        # all triangle edges for if going towards that edge actually takes us
        # closer to p2

        if λp2min ≥ -interpolate_eps
            SI.integrator(summand, ((p1, p2)), (bp1, bp2), icell) #!!! allocates; probably due to .integrator being Any
            result += summand
            break
        else
            # calculate intersection point with corresponding edge
            imin = 0
            t = 1.0 + interpolate_eps
            pint .= p1
            closestimin = 1
            closestdist = Inf
            p1_temp3 .= 0

            while !iscloser(pint, p1, p2, interpolate_eps) ||
                    (any(x -> x <= -interpolate_eps, bpint) || any(x -> x >= 1 + interpolate_eps, bpint)) ||
                    dist(bpint, bp1) <= interpolate_eps #!!! .< allocates a temp vector
                # check if pint takes us closer to p2 and if bpint is inside the cell *and* if we actually moved with pint from p1
                imin += 1
                if imin == 4
                    imin = closestimin
                    bpint .= p1_temp3
                    break
                end
                t = bp1[imin] / (bp1[imin] - bp2[imin])
                bpint = bp1 + t * (bp2 - bp1)
                eval_trafo!(pint, L2G, bpint)
                if dist(pint, p2) < closestdist && ((all(x -> x >= -interpolate_eps, bpint) && all(x -> x <= 1 + interpolate_eps, bpint)))
                    closestimin = imin
                    closestdist = dist(pint, p2)
                    p1_temp3 .= bpint
                end
            end
            eval_trafo!(pint, L2G, bpint)
            SI.integrator(summand, ((p1, pint)), (bp1, bpint), icell)
            result += summand

            # proceed to next cell along edge of smallest barycentric coord
            prevcell = icell
            icell = xFaceCells[1, xCellFaces[facetogo[1][imin], icell]] #!!! allocates
            icell = icell == prevcell ? xFaceCells[2, xCellFaces[facetogo[1][imin], icell]] : icell #!!! allocates

            p1 .= pint
        end
    end

    if axisymmetric
        return dot(result, hnormal) / (avg_r * edge_length)
    else
        return dot(result, hnormal) / edge_length
    end
end

end
