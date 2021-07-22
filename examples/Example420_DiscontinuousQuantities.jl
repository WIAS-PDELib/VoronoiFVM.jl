#=

# 420: Discontinuous Quantities
 ([source code](SOURCE_URL))

Test  jumping species and quantity handling

=#
module Example420_DiscontinuousQuantities
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids
using GridVisualize
using LinearAlgebra


function main(;N=5, Plotter=nothing,unknown_storage=:sparse)
    XX=collect(0:0.1:1)
    xcoord=XX
    for i=1:N-1
	xcoord=glue(xcoord,XX.+i)
    end
    grid2=simplexgrid(xcoord)
    for i=1:N
	cellmask!(grid2,[i-1],[i],i)
    end	
    for i=1:N-1
	bfacemask!(grid2,[i],[i],i+2)
    end

    params=zeros(2,num_cellregions(grid2))
    for i=1:num_cellregions(grid2)
        params[1,i]=i
        params[2,i]=10*i
    end

    
    system=VoronoiFVM.System(grid2,unknown_storage=unknown_storage)

    ## First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitely.
    cspec=ContinuousQuantity(system,1:N;ispec=1,id=1)

    ## A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user.
    dspec=DiscontinuousQuantity(system,1:N; regionspec=[2+i%2 for i=1:N],id=2)

    carrierList = [cspec dspec]
    numberCarriers = length(carrierList)

    params2=zeros(1, numberCarriers)

    for icc in carrierList
        params2[icc] = 2
    end

    for i=1:num_cellregions(grid2)
        @assert params[cspec,i] == i
        @assert params[dspec,i] == 10*i
    end

    for i=1:num_cellregions(grid2)
        params[cspec,i] = -i
        params[dspec,i] = -10*i
    end

    for i=1:num_cellregions(grid2)
        @assert params[1,i] == -i
        @assert params[2,i] == -10*i
    end

    
    ##For both quantities, we define simple diffusion fluxes:

    function flux(f,u,edge)
	f[dspec]=u[dspec,1]-u[dspec,2]
	f[cspec]=u[cspec,1]-u[cspec,2]
    end


    d1=1
    q1=0.2

    ## Define a thin layer inteface condition for `dspec` and an interface source for `cspec`.
    function breaction(f,u,node)
	if node.region>2
            react=(u[dspec,1]-u[dspec,2])/d1
            f[dspec,1]= react
            f[dspec,2]= -react
	    f[cspec]=-q1*u[cspec]
        end
    end
    

    physics!(system,VoronoiFVM.Physics(
        flux=flux,
        breaction=breaction
    ))
    
    
    ## Set boundary conditions
    boundary_dirichlet!(system,dspec,1,1)
    boundary_dirichlet!(system,dspec,2,0.1)
    boundary_dirichlet!(system,cspec,1,0.1)
    boundary_dirichlet!(system,cspec,2,1.0)
    subgrids=VoronoiFVM.subgrids(dspec,system)


    U=solve(unknowns(system,inival=0),system)

    dvws=VoronoiFVM.views(U,dspec,subgrids,system)
    cvws=VoronoiFVM.views(U,cspec,subgrids,system)
    vis=GridVisualizer(resolution=(600,300),Plotter=Plotter)
    for i=1:length(dvws)
	scalarplot!(vis,subgrids[i],dvws[i],flimits=(-0.5,1.5),clear=false,color=:red)
	scalarplot!(vis,subgrids[i],cvws[i],flimits=(-0.5,1.5),clear=false,color=:green)
    end
    reveal(vis)
    norm(system,U,2)
end


function test()
    testval=8.884764592190459
    main(unknown_storage=:sparse) ≈ testval &&
        main(unknown_storage=:dense) ≈ testval
end

end
