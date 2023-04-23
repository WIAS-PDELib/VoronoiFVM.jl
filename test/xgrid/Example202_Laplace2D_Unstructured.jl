#=

# 202: 2D Laplace equation on unstructured mesh 
([source code](SOURCE_URL))

=#
module Example202_Laplace2D_Unstructured

using VoronoiFVM
using LinearAlgebra


# Flux function which describes the flux
# between neighboring control volumes $\omega_k$ and $\omega_l$
function flux!(f,u,edge)
    f[1]=u[1,1]-u[1,2]
end

function source!(f,edge)
    f[1]=1
end


function main(;Plotter=nothing, plot_grid=false,nref=0)

    nspecies=1
    ispec=1
    nrad=10*2^nref
    grid=VoronoiFVM.Grid(points=reduce(hcat,[ [cos(2*pi*i/nrad), sin(2*pi*i/nrad)] for i=1:nrad]),
                     bfaces=reduce(hcat,vcat([ [i,i+1] for i=1:nrad-1],[[nrad,1]])),                  
                     bfaceregions=ones(nrad),
                     regionpoints=[0.0 0.0;],
                     regionnumbers=[1],
                     regionvolumes=[0.1*2.0^(-2*nref)])

    if plot_grid
        p=plot(Plotter=Plotter,grid)
        return
    end

    
    # Create a physics structure
    physics=VoronoiFVM.Physics(num_species=nspecies,flux=flux!,source=source!)

    # Create a finite volume system with dense storage of unknowns
    sys=VoronoiFVM.System(grid,physics,unknown_storage=:dense)

    # Enable species 1 in region 1
    enable_species!(sys,ispec,[1])

    # Set boundary conditions
    # Dirichlet boundary conditions are marked by setting a corresponding value of the
    # boundary factor
    for i=1:num_bfaceregions(grid)
        boundary_dirichlet!(sys,ispec,i,0.0)
    end
    
    # Create & initialize array for solution and initial value
    inival=unknowns(sys)
    solution=unknowns(sys)
    inival.=0
    solution.=0

    # Solve stationary problem
    solve!(solution,inival,sys)

    if Plotter!=nothing
        plot(grid, solution[1,:], Plotter=Plotter)
    end
    
    # Return test value
    return norm(solution,Inf)
end

#
# Called by unit test
#
function test()
    main() ≈ 0.23422190956113936
end

end
