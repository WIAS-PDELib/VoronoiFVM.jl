# # 120: Differing species sets in regions, 1D
# ([source code](SOURCE_URL))

module Example120_ThreeRegions1D

using Printf
using VoronoiFVM




function main(;n=30,Plotter=nothing,plot_grid=false, verbose=false,dense=false)
    h=3.0/(n-1)
    X=collect(0:h:3.0)
    grid=VoronoiFVM.Grid(X)
    cellmask!(grid,[0.0],[1.0],1)
    cellmask!(grid,[1.0],[2.1],2)
    cellmask!(grid,[1.9],[3.0],3)

    subgrid1=subgrid(grid,[1])   
    subgrid2=subgrid(grid,[1,2,3])
    subgrid3=subgrid(grid,[3])
    
    if isplots(Plotter)&&plot_grid
        p=Plotter.plot()
        VoronoiFVM.plot(Plotter,grid,p=p)
        Plotter.gui(p)
        return
    end
    
    eps=[1,1,1]
    k=[1,1,1]

    physics=VoronoiFVM.Physics(
    num_species=3,
    reaction=function(f,u,node,data)
        if node.region==1
            f[1]=k[1]*u[1]
            f[2]=-k[1]*u[1]
        elseif node.region==3
            f[2]=k[3]*u[2]
            f[3]=-k[3]*u[2]
        else
            f[1]=0
        end
    end,
    
    flux=function(f,u,edge,data)   
        uk=viewK(edge,u)
        ul=viewL(edge,u)
        if edge.region==1
            f[1]=eps[1]*(uk[1]-ul[1])
            f[2]=eps[2]*(uk[2]-ul[2])
        elseif edge.region==2
            f[2]=eps[2]*(uk[2]-ul[2])
        elseif edge.region==3
            f[2]=eps[2]*(uk[2]-ul[2])
            f[3]=eps[3]*(uk[3]-ul[3])
        end
    end,
    
    source=function(f,node,data)
        if node.region==1
            f[1]=1.0e-4*(3.0-node.coord[1])
        end
    end,
    
    storage=function(f,u,node,data)
        f.=u
    end
    )

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end

    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1,2,3])
    enable_species!(sys,3,[3])

    boundary_dirichlet!(sys,3,2,0.0)
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    tstep=0.01
    time=0.0
    istep=0
    tend=10
    testval=0
    while time<tend
        time=time+tstep
        solve!(U,inival,sys,control=control,tstep=tstep)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        tstep*=1.0
        istep=istep+1
        testval=U[2,15]
        if isplots(Plotter)
            Plots=Plotter
            p=Plots.plot()
            VoronoiFVM.plot(Plots,subgrid1, U[1,:],label="spec1", color=(0.5,0,0),p=p,show=false)
            VoronoiFVM.plot(Plots,subgrid2, U[2,:],label="spec2", color=(0.0,0.5,0),p=p,show=false)
            VoronoiFVM.plot(Plots,subgrid3, U[3,:],label="spec3", color=(0.0,0.0,0.5),p=p,show=false)
            Plots.gui(p)
        end
    end
    return testval
end

function test()
    main(dense=false) ≈ 0.00039500514567080265 &&
        main(dense=true) ≈ 0.00039500514567080265
end

end