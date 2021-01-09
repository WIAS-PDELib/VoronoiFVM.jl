# # 110: 1D Reaction Diffusion equation with two species
# ([source code](SOURCE_URL))
# Solve the nonlinear coupled reaction diffusion problem
# 
# ```math
# -\nabla (0.01+2u_2)\nabla u_1 + u_1u_2= 0.0001(0.01+x)
# ```
# 
# ```math
# -\nabla (0.01+2u_1)\nabla u_2 - u_1u_2 = 0.0001(1.01-x)
# ```
# 
# 
# in $\Omega=(0,1)$ with boundary condition $u_1(0)=1$, $u_2(0)=0$ and $u_1(1)=1$, $u_2(1)=1$.
# 


module Example110_ReactionDiffusion1D_TwoSpecies

using Printf
using VoronoiFVM
using ExtendableGrids

function main(;n=100,Plotter=nothing,verbose=false,unknown_storage=:sparse)
    h=1/n
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
    
    eps=[1.0,1.0]
    
    physics=VoronoiFVM.Physics(num_species=2,
                               
                               reaction=function(f,u,node)
                               f[1]=u[1]*u[2]
                               f[2]=-u[1]*u[2]
                               end,
                               
                               flux=function(f,u0,edge)
                               u=unknowns(edge,u0)

                               nspecies=2
                               f[1]=eps[1]*(u[1,1]-u[1,2])*(0.01+u[2,1]+u[2,2])
                               f[2]=eps[2]*(u[2,1]-u[2,2])*(0.01+u[1,1]+u[1,2])
                               end,
                               
                               source=function(f,node)
                               f[1]=1.0e-4*(0.01+node[1])
                               f[2]=1.0e-4*(0.01+1.0-node[1])
                               end,
                               
                               storage=function(f,u,node)
                               f[1]=u[1]
                               f[2]=u[2]
                               end
                               )
    
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    boundary_dirichlet!(sys,1,1,1.0)
    boundary_dirichlet!(sys,1,2,0.0)
    
    boundary_dirichlet!(sys,2,1,1.0)
    boundary_dirichlet!(sys,2,2,0.0)
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    control.damp_initial=0.1
    u5=0
    p=GridPlotContext(Plotter=Plotter,layout=(2,1))
    for xeps in [1.0,0.5,0.25,0.1,0.05,0.025,0.01]
        eps=[xeps,xeps]
        solve!(U,inival,sys,control=control)
        inival.=U
        gridplot!(p[1,1],grid,U[1,:],clear=true,title="U1, eps=$(xeps)")
        gridplot!(p[2,1],grid,U[2,:],clear=true,title="U2, eps=$(xeps)",reveal=true)
        sleep(0.2)
        u5=U[5]
    end
    return u5
end

function test()
    testval=0.7117546972922056
    main(unknown_storage=:sparse) ≈ testval && main(unknown_storage=:dense) ≈ testval
end
end

