#=

# 107: 1D Nonlinear Storage 
([source code](SOURCE_URL))

This equation comes from the transformation of the nonlinear diffuision equation.
```math
\partial_t u^\frac{1}{m} -\Delta u = 0
```
in $\Omega=-1,1)$ with homogeneous Neumann boundary conditons.
We can derive an exact solution from the Barenblatt solution of the previous
example.

=#

module Example107_NonlinearStorage1D
using Printf
using VoronoiFVM
using ExtendableGrids

function barenblatt(x,t,m)
    tx=t^(-1.0/(m+1.0))
    xx=x*tx
    xx=xx*xx
    xx=1- xx*(m-1)/(2.0*m*(m+1));
    if xx<0.0
        xx=0.0
    end
    return tx*xx^(1.0/(m-1.0))
end


function main(;n=20,m=2.0,Plotter=nothing,verbose=false, unknown_storage=:sparse,tend=0.01,tstep=0.0001)
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n/2)
    X=collect(-1:h:1)
    grid=VoronoiFVM.Grid(X)

    ## Flux function which describes the flux
    ## between neigboring control volumes
    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=u[1,1]-u[1,2]
    end

    ϵ=1.0e-10

    ## Storage term
    ## This needs to be regularized as its derivative
    ## at 0 is infinity
    function storage!(f,u,node)
        f[1]=(ϵ+u[1])^(1.0/m)
    end
    
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
    
    
    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    ## Add species 1 to region 1
    enable_species!(sys,1,[1])
    
    ## Create a solution array
    inival=unknowns(sys)
    solution=unknowns(sys)
    t0=0.001

    ## Broadcast the initial value
    inival[1,:].=map(x->barenblatt(x,t0,m)^m,X)
    solution.=inival

    ## Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    time=t0
    p=GridPlotContext(Plotter=Plotter,layout=(2,1))
    while time<tend
        time=time+tstep
        solve!(solution,inival,sys,control=control,tstep=tstep)
        inival.=solution
        gridplot!(p[1,1],grid,solution[1,:],title=@sprintf("numerical, t=%.5f",time),clear=true)
        gridplot!(p[2,1],grid,map(x->barenblatt(x,time,m)^m,grid),title=@sprintf("exact, t=%.4f",time),clear=true)
        sleep(1.0e-5)
#        yield()
    end
    return sum(solution)
end


function test()
    testval=173.84998139534142
    main(unknown_storage=:sparse) ≈ testval && main(unknown_storage=:dense) ≈ testval
end

end 

