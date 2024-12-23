module test_geomitems
using VoronoiFVM
using ExtendableGrids
using ForwardDiff
using Test

const Dual64 = ForwardDiff.Dual{Float64, Float64, 1}

flux(y, u, edge, data) = y[1] = u[1, 1] - u[1, 2]


function main(;n=3)
    g = simplexgrid(0:0.1:1)
    sys = VoronoiFVM.System(g; flux,species = 1:n)
    time=0.0
    λ=0.0
    params=zeros(3)

    utest=zeros(Dual64,2*n)
    node = VoronoiFVM.Node(sys, time, λ, params)
    u=unknowns(node,utest)
    @test typeof(u[1])==eltype(u)
    @test typeof(copy(u)[1])==eltype(u)
    @test length(u)==n
    @test size(u)==(n,)

    bnode = VoronoiFVM.BNode(sys, time, λ, params)
    u=unknowns(bnode,utest)
    @test typeof(u[1])==eltype(u)
    @test typeof(copy(u)[1])==eltype(u)
    @test length(u)==n
    @test size(u)==(n,)

    edge = VoronoiFVM.Edge(sys, time, λ, params)
    u=unknowns(edge,utest)
    @test typeof(u[1])==eltype(u)
    @test typeof(copy(u)[1])==eltype(u)
    @test length(u)==2*n
    @test size(u)==(n,2)

    bedge = VoronoiFVM.BEdge(sys, time, λ, params)
    u=unknowns(bedge,utest)
    @test typeof(u[1])==eltype(u)
    @test typeof(copy(u)[1])==eltype(u)
    @test length(u)==2*n
    @test size(u)==(n,2)



    ftest=zeros(Dual64,2*n)
    node = VoronoiFVM.Node(sys, time, λ, params)
    f=VoronoiFVM.rhs(node,ftest)
    @test typeof(f[1])==eltype(f)
    @test typeof(copy(f)[1])==eltype(f)
    @test length(f)==n
    @test size(f)==(n,)

    bnode = VoronoiFVM.BNode(sys, time, λ, params)
    f=VoronoiFVM.rhs(bnode,ftest)
    @test typeof(f[1])==eltype(f)
    @test typeof(copy(f)[1])==eltype(f)
    @test length(f)==n
    @test size(f)==(n,)

    edge = VoronoiFVM.Edge(sys, time, λ, params)
    f=VoronoiFVM.rhs(edge,ftest)
    @test typeof(f[1])==eltype(f)
    @test typeof(copy(f)[1])==eltype(f)
    @test length(f)==n
    @test size(f)==(n,)

    bedge = VoronoiFVM.BEdge(sys, time, λ, params)
    f=VoronoiFVM.rhs(bedge,ftest)
    @test typeof(f[1])==eltype(f)
    @test typeof(copy(f)[1])==eltype(f)
    @test length(f)==n
    @test size(f)==(n,)


end

function runtests()
    main()
    return nothing
end
end
