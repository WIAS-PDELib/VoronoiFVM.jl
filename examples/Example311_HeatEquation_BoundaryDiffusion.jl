#=

# 311: Heat Equation with boundary diffusion 
([source code](SOURCE_URL))


=#

module Example311_HeatEquation_BoundaryDiffusion
using Printf
using VoronoiFVM
using ExtendableGrids

"""
  We solve the following system

      ∂_tu - εΔu = 0            in [0,T] × Ω
           ε∇u⋅ν = k(u-v)       on [0,T] × Γ_1
           ε∇u⋅ν = 0            on [0,T] × (∂Ω ∖ Γ_1)
  ∂_tv -ε_ΓΔ_Γ v = f(x) +k(u-v) on [0,T] × Γ_1
          u(0)   = 0.5          in   {0} × Ω
          v(0)   = 0.5          on   {0} × Γ_1  
"""

function main(n=1)

    breg = 5 # boundary region number for surface diffusion

    println("Surface diffusion example")

    hmin   = 0.01 * 2.0^(-n+1)
    hmax   = 0.1  * 2.0^(-n+1)
    XLeft  = geomspace(0.0, 0.5, hmax, hmin)
    XRight = geomspace(0.5, 1.0, hmin, hmax)
    X      = glue(XLeft, XRight)
    
    Z      = geomspace(0.0, 1.0, hmin, 2*hmax)
    
    grid   = VoronoiFVM.Grid(X, X, Z)


    # parameters
    eps      = 1.0e0  # bulk heat conduction coefficient
    eps_surf = 1.0e-2 # surface diffusion coefficient
    k        = 1.0    # transmission coefficient
    physics  = VoronoiFVM.Physics(
        num_species = 2,
        flux        = function(f,u0,edge)
                            u    = unknowns(edge,u0)
                            f[1] = eps*(u[1,1]-u[1,2])
                      end,
        
        bflux       = function(f, u0, edge)
                        if edge.region == breg
                          u        = unknowns(edge,u0)
                          f[2]     = eps_surf*(u[2,1]-u[2,2])
                        else
                          f[2] = 0.0
                        end

                      end,
        breaction   = function(f,u,node)
                        if node.region == breg
                          f[1] = k*(u[1]-u[2])
                          f[2] = k*(u[2]-u[1])
                        else
                          f[1] = 0.0        
                          f[2] = 0.0
                        end
                      end,
        bsource      = function(f, bnode)
                        x1   = bnode[1]-0.5
                        x2   = bnode[2]-0.5
                        x3   = bnode[3]-0.5
                        f[2] = 1.0e4*exp(-20.0*(x1^2+x2^2+x3^2))
                      end,
 
        bstorage    = function(f,u,node)
                        if node.region == breg
                          f[2] = u[2]
                        end
                      end,

        storage     = function(f,u,node)
                        f[1] = u[1]
                      end
    )

    sys=VoronoiFVM.System(grid,physics,unknown_storage=:sparse)
    enable_species!(sys,1,[1])
    enable_boundary_species!(sys,2,[breg])

    function tran32!(a,b)
        a[1]=b[2]
    end
    
    bgrid2=subgrid(grid,[breg],boundary=true,transform=tran32!)
    
    inival  = unknowns(sys)
    U       = unknowns(sys)
    inival .= 0.5


    control             = VoronoiFVM.NewtonControl()
    control.verbose     = true
    control.tol_linear  = 1.0e-5
    control.max_lureuse = 10

    tstep = 0.1
    time  = 0.0
    step  = 0
    T     = 1.0
    while time < T
        time=time+tstep
        solve!(U, inival, sys, control=control, tstep=tstep)
        inival.=U
        U_surf=view(U[2,:],bgrid2)

        if true
            @printf("time=%g\n",time)
        end

        tstep*=1.0

        step += 1
    end

    nothing
end


function test()
    main()
    true
end

end