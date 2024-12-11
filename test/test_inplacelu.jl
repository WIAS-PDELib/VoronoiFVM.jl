module test_inplacelu
using VoronoiFVM
using Test
using StrideArraysCore: StrideArray, @gc_preserve, StaticInt
using StaticArrays
using Random, LinearAlgebra
using ForwardDiff

"""
    Test inplace linear system solution via non-pivoting inplace_linsolve!
    which implements Doolittle's method in VoronoiFVM. 
    Uses StrideArraysCore.StrideArray.
    This must not allocate.
"""
function inplacelu_nopiv_stridearray(::Type{Val{N}}, ::Type{T}) where {N, T}
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

inplacelu_nopiv_stridearray(n, T) = inplacelu_nopiv_stridearray(Val{n}, T)

"""
    Test inplace linear system solution via non-pivoting inplace_linsolve!
    which implements Doolittle's method in VoronoiFVM. 
    Uses StaticArrays.MArray.
    This must not allocate.
"""
function inplacelu_nopiv_marray(::Type{Val{N}}, ::Type{T}) where {N, T}
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

inplacelu_nopiv_marray(n, T) = inplacelu_nopiv_marray(Val{n}, T)


"""
    Test inplace linear system solution via pivoting inplace_linsolve!
    which is implemented using RecursiveFactorization.lu
    Uses StrideArraysCore.StrideArray.
    This must not allocate.
"""
function inplacelu_piv_stridearray(::Type{Val{N}}, ::Type{T}) where {N, T}
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

inplacelu_piv_stridearray(n, T) = inplacelu_piv_stridearray(Val{n}, T)

"""
    Test inplace linear system solution via pivoting inplace_linsolve!
    which is implemented using RecursiveFactorization.lu
    Uses StaticArrays.MArray.
    This must not allocate.
"""
function inplacelu_piv_marray(::Type{Val{N}}, ::Type{T}) where {N, T}
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

inplacelu_piv_marray(n, T) = inplacelu_piv_marray(Val{n}, T)

# Define dual number type
const Dual64 = ForwardDiff.Dual{Float64, Float64, 1}


function runtests()
    # Check if Dual64 is fully parametrized
    @test isbitstype(Dual64)
    
    # first precompile to avoid allocations during precompilation
    inplacelu_nopiv_marray(10, Float64)
    inplacelu_nopiv_stridearray(10, Float64)
    inplacelu_nopiv_marray(10, Dual64)
    inplacelu_nopiv_stridearray(10, Dual64)

    n1 = @allocated inplacelu_nopiv_stridearray(10, Float64)
    @test n1 == 0
    n2 = @allocated inplacelu_nopiv_marray(10, Float64)
    @test n2 == 0
    n3 = @allocated inplacelu_nopiv_marray(10, Dual64)
    @test n3 == 0
    n4 = @allocated inplacelu_nopiv_stridearray(10, Dual64)
    @test n4 == 0

    inplacelu_piv_stridearray(10, Float64)
    inplacelu_piv_marray(10, Float64)
    inplacelu_piv_stridearray(10, Dual64)
    inplacelu_piv_marray(10, Dual64)

    m1 = @allocated inplacelu_piv_stridearray(10, Float64)
    @test m1 == 0
    m2 = @allocated inplacelu_piv_marray(10, Float64)
    @test m2 == 0
    m3 = @allocated inplacelu_piv_marray(10, Dual64)
    @test m3 == 0
    m4 = @allocated inplacelu_piv_stridearray(10, Dual64)
    @test m4 == 0
    return true
end

end
