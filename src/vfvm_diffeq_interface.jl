#
# Interface to VoronoiFVMDiffEq.jl
#
"""
    $(TYPEDEF)

History information for DiffEqHistory
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct DiffEqHistory 
    """ number of combined jacobi/rhs evaluations"""
    nd=0
    """ number of combined jacobi evaluations"""
    njac=0
    """ number of rhs evaluations"""
    nf=0
end

details(h::DiffEqHistory)=(nd=h.nd,njac=h.njac,nf=h.nf)
Base.summary(h::DiffEqHistory)=details(h)


#Evaluate functiaon and Jacobian at u if they have not been evaluated before at u.
#See https://github.com/SciML/DifferentialEquations.jl/issues/521 for discussion of another way to do this.
function _eval_res_jac!(sys,u,t)
    uhash=hash(u)
    if uhash!=sys.uhash
        ur=reshape(u,sys)
        eval_and_assemble(sys,ur,ur,sys.residual,value(t),Inf,0.0,zeros(0))
        sys.uhash=uhash
        sys.history.nd+=1
    end
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
rhs function for DifferentialEquations.jl.
"""
function eval_rhs!(du, u, sys,t)
    _eval_res_jac!(sys,u,t)
    du.=-vec(sys.residual)
    sys.history.nf+=1
    nothing
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
jacobi matrix calculation function for DifferentialEquations.jl.
"""
function eval_jacobian!(J, u, sys,t)
    _eval_res_jac!(sys,u,t)
    # Need to implement broadcast for ExtendableSparse.
    J.=-sys.matrix.cscmatrix
    sys.history.njac+=1
    nothing
end

"""
$(SIGNATURES)

Calculate the mass matrix for use with DifferentialEquations.jl.
Return a Diagonal matrix if it occurs to be diagonal, otherwise return a SparseMatrixCSC.
"""
function mass_matrix(system::AbstractSystem{Tv, Tc, Ti, Tm}) where {Tv, Tc, Ti, Tm}
    physics=system.physics
    grid=system.grid
    data=physics.data
    node=Node(system)
    bnode=BNode(system)
    nspecies=num_species(system)
    cellnodefactors=system.cellnodefactors
    bfacenodefactors=system.bfacenodefactors
    bgeom=grid[BFaceGeometries][1]
    nbfaces=num_bfaces(grid)
    ndof=num_dof(system)

    UK=Array{Tv,1}(undef,nspecies)
    stor_eval = ResJacEvaluator(physics,:storage,UK,node,nspecies)
    bstor_eval = ResJacEvaluator(physics,:bstorage,UK,node,nspecies)

    U=unknowns(system,inival=0)
    M=ExtendableSparseMatrix{Tv,Tm}(ndof,ndof)
    
    geom=grid[CellGeometries][1]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    bfacenodes=grid[BFaceNodes]

    for icell=1:num_cells(grid)
        ireg=cellregions[icell]
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)
            K=node.index
            for ispec=1:nspecies
                UK[ispec]=U[ispec,K]
            end

            evaluate!(stor_eval,UK)
            res_stor=res(stor_eval)
            jac_stor=jac(stor_eval)

            for idof=_firstnodedof(U,K):_lastnodedof(U,K)
                ispec=_spec(U,idof,K)
                if system.region_species[ispec,ireg]<=0
                    continue
                end
                for jdof=_firstnodedof(U,K):_lastnodedof(U,K)
                    jspec=_spec(U,jdof,K)
                    if system.region_species[jspec,ireg]<=0
                        continue
                    end
                    _addnz(M,idof,jdof,jac_stor[ispec,jspec],cellnodefactors[inode,icell])
                end
            end
        end
    end
    
    if docall(bstor_eval)
        for ibface=1:nbfaces
            for ibnode=1:num_nodes(bgeom)
                _fill!(bnode,ibnode,ibface)
                K=bnode.index
                for ispec=1:nspecies
                    UK[ispec]=U[ispec,bnode.index]
                end
                evaluate!(bstor_eval,UK)
                res_bstor=res(bstor_eval)
                jac_bstor=jac(bstor_eval)
                for idof=_firstnodedof(U,K):_lastnodedof(U,K)
                    ispec=_spec(U,idof,K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    for jdof=_firstnodedof(U,K):_lastnodedof(U,K)
                        jspec=_spec(U,jdof,K)
                        if !isdof(system,jspec,K)
                            continue
                        end
                        _addnz(M,idof,jdof,jac_bstor[ispec,jspec],bfacenodefactors[ibnode,ibface])
                    end
                end
                
            end
        end
    end
    Mcsc=sparse(M)
    isdiag(Mcsc) ? Diagonal([Mcsc[i,i] for i=1:ndof]) : Mcsc
end

"""
$(SIGNATURES)
Prepare system for use with VoronoiFVMDiffEq.

- `jacval`: value at which to evaluate jacobian to obtatin prototype
- `tjac`: time moment for jacobian
 
Returns a prototype for the jacobian.
"""
function prepare_diffeq!(sys,jacval,tjac)
    sys.history=DiffEqHistory()
    _complete!(sys, create_newtonvectors=true)
    _eval_res_jac!(sys,jacval,tjac)
    flush!(sys.matrix)
    sys.matrix.cscmatrix
end
