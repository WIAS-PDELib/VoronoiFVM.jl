# # 220: 2D Nonlinear Poisson with boundary reaction and boundary species
# ([source code](SOURCE_URL))

module Example220_NonlinearPoisson2D_BoundarySpecies

using Printf
using VoronoiFVM
using ExtendableGrids

function main(;n=10,Plotter=nothing,verbose=false,unknown_storage=:sparse)
    
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    


    grid=VoronoiFVM.Grid(X,Y)
    
    
    k=1.0
    eps=1.0
    physics=VoronoiFVM.Physics(
    num_species=3,
    breaction=function(f,u,node)
        if  node.region==2
            f[1]=k*(u[1]-u[3])
            f[3]=k*(u[3]-u[1])+ k*(u[3]-u[2])
            f[2]=k*(u[2]-u[3])
        end
    end,
    
    bstorage=function(f,u,node)
        if  node.region==2
            f[3]=u[3]
        end
    end,
    
    
    flux=function(f,u0,edge)
        u=unknowns(edge,u0)
        f[1]=eps*(u[1,1]-u[1,2])
        f[2]=eps*(u[2,1]-u[2,2])
    end,
    
    source=function(f,node)
        x1=node[1]-0.5
        x2=node[2]-0.5
        f[1]=exp(-20.0*(x1^2+x2^2))
    end,
    
    storage=function(f,u,node)
        f[1]=u[1]
        f[2]=u[2]
    end
    )
    
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])
    enable_boundary_species!(sys,3,[2])

    
    function tran32!(a,b)
        a[1]=b[2]
    end
    
    bgrid2=subgrid(grid,[2],boundary=true,transform=tran32!)
   
    inival=unknowns(sys)
    inival.=0.0
    U=unknowns(sys)

    eps=1.0e-2
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.tol_linear=1.0e-5
    control.tol_relative=1.0e-5
    control.max_lureuse=0
    tstep=0.01
    time=0.0
    istep=0
    u5=0
    p=GridPlotContext(Plotter=Plotter,layout=(3,1))
    while time<1
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        U_bound=view(U[3,:],bgrid2)
        u5=U_bound[5]
        gridplot!(p[1,1],grid,U[1,:],clear=true)
        gridplot!(p[2,1],grid,U[2,:])
        gridplot!(p[3,1],bgrid2,U_bound,show=true,flimits=(0,0.0025))
    end
    return u5
end

function test()
    main(unknown_storage=:sparse) ≈ 0.0020781361856598
    main(unknown_storage=:dense) ≈ 0.0020781361856598
end
end
