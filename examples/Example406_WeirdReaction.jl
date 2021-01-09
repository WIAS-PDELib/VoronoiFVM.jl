#=

# 406: 1D Weird Surface Reaction
 ([source code](SOURCE_URL))


Species $A$ and $B$ exist in the interior of the domain, species $C$
lives a the boundary $\Gamma_1$.  We assume a heterogeneous reaction scheme
where $A$ reacts to $B$ with a rate depending on $\nabla A$ near the surface

```math
\begin{aligned}
      A &\leftrightarrow B\\
\end{aligned}
```

In $\Omega$, both $A$ and $B$ are transported through diffusion:

```math
\begin{aligned}
\partial_t u_B - \nabla\cdot D_A \nabla u_A & = f_A\\
\partial_t u_B - \nabla\cdot D_B \nabla u_B & = 0\\
\end{aligned}
```
Here, $f(x)$ is a source term creating $A$.
On $\Gamma_2$, we set boundary conditions
```math
\begin{aligned}
D_A \nabla u_A & = 0\\
u_B&=0
\end{aligned}
```
describing no normal flux for $A$ and zero concentration of $B$.
On $\Gamma_1$, we use the mass action law to describe the boundary reaction and
the evolution of the boundary concentration $C$. We assume that there is a limited
amount of surface sites $S$ for species C, so in fact A has to react with a free
surface site in order to become $C$ which reflected by the factor $1-u_C$. The same
is true for $B$.
```math
\begin{aligned}
R_{AB}(u_A, u_B)&=k_{AB}^+exp(u_A'(0))u_A - k_{AB}^-exp(-u_A'(0))u_B\\
- D_A \nabla u_A  +  R_{AB}(u_A, u_B)& =0 \\
- D_B \nabla u_B  -  R_{AB}(u_A, u_B)& =0 \\
\end{aligned}
```

=#
module Example406_WeirdReaction
using Printf
using VoronoiFVM
using SparseArrays
using ExtendableGrids

function main(;n=10,
              Plotter=nothing,
              verbose=false,
              tend=1,
              unknown_storage=:sparse,
              autodetect_sparsity=true)
    
    h=1.0/convert(Float64,n)
    X=collect(0.0:h:1.0)
    N=length(X)
    
    grid=VoronoiFVM.Grid(X)
    ## By default, \Gamma_1 at X[1] and \Gamma_2 is at X[end]
    
    ## Species numbers
    iA=1
    iB=2
    iC=3


    ## Diffusion flux for species A and B
    D_A=1.0
    D_B=1.0e-2
    function flux!(f,u0,edge)
        u=unknowns(edge,u0)
        f[iA]=D_A*(u[iA,1]-u[iA,2])
        f[iB]=D_B*(u[iB,1]-u[iB,2])
    end

    ## Storage term of species A and B
    function storage!(f,u,node)
        f[iA]=u[iA]
        f[iB]=u[iB]
    end

    ## Source term for species a around 0.5
    function source!(f,node)
        x1=node[1]-0.5
        f[iA]=exp(-100*x1^2)
    end

    ## Reaction constants (p = + , m = -)
    ## Choosen to prefer path A-> B
    kp_AB=1.0
    km_AB=0.1
    
    
    function breaction!(f,u,node)
        if  node.region==1
            R=kp_AB*exp(u[iC])*u[iA] - exp(-u[iC])*km_AB*u[iB]
            f[iA]+=R
            f[iB]-=R
        end
    end

    ## This generic operator works on the full solution seen as linear vector, and indexing
    ## into it needs to be performed with the help of idx (defined below for a solution vector)
    ## Its sparsity is detected automatically using SparsityDetection.jl
    ## Here, we calculate the gradient of u_A at the boundary and store the value in u_C which
    ## is then used as a parameter in the boundary reaction
    function generic_operator!(f,u,sys)
        f.=0
        f[idx[iC,1]]=u[idx[iC,1]]  + 0.1*(u[idx[iA,1]]-u[idx[iA,2]])/(X[2]-X[1])
    end

    # If we know the sparsity pattern, we can here create a
    # sparse matrix with values set to 1 in the nonzero
    # slots. This allows to circumvent the
    # autodetection which may takes some time.
    function generic_operator_sparsity(sys)
        idx=unknown_indices(unknowns(sys))
        sparsity=spzeros(num_dof(sys),num_dof(sys))
        sparsity[idx[iC,1],idx[iC,1]]=1
        sparsity[idx[iC,1],idx[iA,1]]=1
        sparsity[idx[iC,1],idx[iA,2]]=1
        sparsity
    end



    
    if autodetect_sparsity
        physics=VoronoiFVM.Physics(
            num_species=3,
            breaction=breaction!,
            generic=generic_operator!,
            flux=flux!,
            storage=storage!,
            source=source!
        )
    else
        physics=VoronoiFVM.Physics(
            num_species=3,
            breaction=breaction!,
            generic=generic_operator!,
            generic_sparsity=generic_operator_sparsity,
            flux=flux!,
            storage=storage!,
            source=source!
        )
    end
    
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    ## Enable species in bulk resp
    enable_species!(sys,iA,[1])
    enable_species!(sys,iB,[1])

    ## Enable surface species
    enable_boundary_species!(sys,iC,[1])
    
    ## Set Dirichlet bc for species B on \Gamma_2
    boundary_dirichlet!(sys,iB,2,0.0)

    ## Initial values
    inival=unknowns(sys)
    inival.=0.0
    U=unknowns(sys)
    idx=unknown_indices(U)
    
    tstep=0.01
    time=0.0
    T=Float64[]
    u_C=Float64[]
    
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    p=GridPlotContext(Plotter=Plotter,layout=(2,1))
    while time<tend
        time=time+tstep
        solve!(U,inival,sys,tstep=tstep,control=control)
        inival.=U
        if verbose
            @printf("time=%g\n",time)
        end
        ## Record  boundary pecies
        push!(T,time)
        push!(u_C,U[iC,1])

        gridplot!(p[1,1],grid,U[iA,:],label="[A]",title=@sprintf("max_A=%.5f max_B=%.5f u_C=%.5f",maximum(U[iA,:]),maximum(U[iB,:]),u_C[end]),color=:red)
        gridplot!(p[1,1],grid,U[iB,:], label="[B]",clear=false,color=:blue)
        gridplot!(p[2,1],copy(T),copy(u_C),label="[C]",clear=true,show=true)
    end
    return U[iC,1]
end

function test()
    testval=0.007027597470502758
    main(unknown_storage=:sparse) ≈ testval &&
        main(unknown_storage=:dense) ≈ testval &&
        main(unknown_storage=:sparse,autodetect_sparsity=false) ≈ testval &&
        main(unknown_storage=:dense,autodetect_sparsity=false) ≈ testval
end

end
