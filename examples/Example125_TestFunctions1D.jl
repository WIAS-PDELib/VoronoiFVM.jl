# # 125: Terminal flux calculation via test functions
# ([source code](SOURCE_URL))

module Example125_TestFunctions1D 
using Printf
using VoronoiFVM
using ExtendableGrids


function main(;n=100,Plotter=nothing,verbose=false,unknown_storage=:sparse)
    h=1/n
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
        
    eps=[1,1.0e-1]

    physics=VoronoiFVM.Physics(
         num_species=2,
    reaction=function(f,u,node)
        f[1]=10*(u[1]-u[2])
        f[2]=10*(u[2]-u[1])
    end,

    flux=function(f,u0,edge)   
        u=unknowns(edge,u0)
        f[1]=eps[1]*(u[1,1]-u[1,2])
        f[2]=eps[2]*(u[2,1]-u[2,2])
    end,
    
    
    storage=function(f,u,node)
        f[1]=u[1]
        f[2]=u[2]
    end
    )
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    boundary_neumann!(sys,1,1,0.01)
    boundary_dirichlet!(sys,2,2,0.0)

    factory=VoronoiFVM.TestFunctionFactory(sys)
    tf1=testfunction(factory,[2],[1])
    tf2=testfunction(factory,[1],[2])
    
    
    U=unknowns(sys)
    inival=unknowns(sys)
    inival[2,:].=0.1
    inival[1,:].=0.1
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    I1=0
    p=GridPlotContext(Plotter=Plotter,layout=(2,1))
    for xeps in [1.0,0.1,0.01]
        eps=[xeps,xeps]
        solve!(U,inival,sys,control=control)
        I1=integrate(sys,tf1,U)
        coord=coordinates(grid)
        inival.=U
        gridplot!(p[1,1],grid,U[1,:])
        gridplot!(p[2,1],grid,U[2,:])
        reveal(p)
        u5=U[5]
    end
    return I1[1]
end

function test()
    main(unknown_storage=:sparse) ≈ 0.01 &&
        main(unknown_storage=:dense) ≈ 0.01
end
end
