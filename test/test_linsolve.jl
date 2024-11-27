module test_linsolve
using VoronoiFVM
using Test
using StrideArraysCore: StrideArray, @gc_preserve, StaticInt
using StaticArrays
using Random, LinearAlgebra
using ForwardDiff

const Dual64 = ForwardDiff.Dual{Float64, Float64, 1}

function checklux0(::Type{Val{N}}, ::Type{T}) where {N, T}
    A = StrideArray{T}(undef, StaticInt(N), StaticInt(N))
    x = StrideArray{T}(undef, StaticInt(N))
    b = StrideArray{T}(undef, StaticInt(N))

    for i in 1:N
        for j in 1:N
            A[i, j] = -rand()
        end
        A[i, i] += 100
        x[i] = one(T)
    end

    @gc_preserve mul!(b, A, x)
    @gc_preserve inplace_linsolve!(A, b)

    nm = 0
    for i in 1:N
        nm += (b[i] - x[i])^2
    end
    @assert sqrt(nm) / N < 100.0 * eps(Float64)
    return nothing
end

checklux(n, T) = checklux0(Val{n}, T)

function checklum0(::Type{Val{N}}, ::Type{T}) where {N, T}
    A = MMatrix{N, N, Float64}(undef)
    x = MVector{N, Float64}(undef)
    b = MVector{N, Float64}(undef)

    for i in 1:N
        for j in 1:N
            A[i, j] = -rand()
        end
        A[i, i] += 100
        x[i] = 1
    end

    @gc_preserve mul!(b, A, x)
    @gc_preserve inplace_linsolve!(A, b)

    nm = 0
    for i in 1:N
        nm += (b[i] - x[i])^2
    end

    @assert sqrt(nm) / N < 100.0 * eps(Float64)
    return nothing
end

checklum(n, T) = checklum0(Val{n}, T)


function checklurx0(::Type{Val{N}}, ::Type{T}) where {N, T}
    A = StrideArray{T}(undef, StaticInt(N), StaticInt(N))
    x = StrideArray{T}(undef, StaticInt(N))
    b = StrideArray{T}(undef, StaticInt(N))
    ipiv = StrideArray{Int64}(undef, StaticInt(N))
    for i in 1:N
        for j in 1:N
            A[i, j] = -rand()
        end
        A[i, i] += 100
        x[i] = 1
    end
    @gc_preserve mul!(b, A, x)
    @gc_preserve inplace_linsolve!(A, b, ipiv)

    nm = 0
    for i in 1:N
        nm += (b[i] - x[i])^2
    end

    @assert sqrt(nm) / N < 100.0 * eps(Float64)
    return nothing
end

checklurx(n, T) = checklurx0(Val{n}, T)

function checklurm0(::Type{Val{N}}, ::Type{T}) where {N, T}
    A = MMatrix{N, N, Float64}(undef)
    x = MVector{N, Float64}(undef)
    b = MVector{N, Float64}(undef)
    ipiv = MVector{N, Int64}(undef)
    for i in 1:N
        for j in 1:N
            A[i, j] = -rand()
        end
        A[i, i] += 100
        x[i] = 1
    end
    @gc_preserve mul!(b, A, x)
    @gc_preserve inplace_linsolve!(A, b, ipiv)

    nm = 0
    for i in 1:N
        nm += (b[i] - x[i])^2
    end

    @assert sqrt(nm) / N < 100.0 * eps(Float64)
    return nothing
end

checklurm(n, T) = checklurm0(Val{n}, T)

function runtests()
    checklum(10, Float64)
    checklum(10, Dual64)

    n4 = @allocated checklum(10, Dual64)
    @test n4 == 0

    n2 = @allocated checklum(10, Float64)
    @test n2 == 0

    checklux(10, Float64)
    checklux(10, Dual64)

    n1 = @allocated checklux(10, Float64)
    @test n1 == 0

    n3 = @allocated checklux(10, Dual64)
    @test n3 == 0

    checklurx(10, Float64)
    checklurm(10, Float64)
    checklurx(10, Dual64)
    checklurm(10, Dual64)

    rn1 = @allocated checklurx(10, Float64)
    @test rn1 == 0
    rn2 = @allocated checklurm(10, Float64)
    @test rn2 == 0
    rn3 = @allocated checklurx(10, Dual64)
    @test rn3 == 0
    rn4 = @allocated checklurm(10, Dual64)
    @test rn4 == 0
    return true
end

end
