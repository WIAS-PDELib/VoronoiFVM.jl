##############################################################
const c1 = 1 / 47_900_160
const c2 = -1 / 1_209_600
const c3 = 1 / 30_240
const c4 = -1 / 720
const c5 = 1 / 12
const c6 = -1 / 2

"""
$(SIGNATURES)

Calculation of Bernoulli function via Horner scheme based on Taylor
coefficients around 0.
"""
function bernoulli_horner(x)
    y = x * c1
    y = x * y
    y = x * (c2 + y)
    y = x * y
    y = x * (c3 + y)
    y = x * y
    y = x * (c4 + y)
    y = x * y
    y = x * (c5 + y)
    y = x * (c6 + y)
    return y = 1 + y
end

# Bernoulli thresholds optimized for Float64
const bernoulli_small_threshold = 0.25
const bernoulli_large_threshold = 50.0
##############################################################
"""
$(SIGNATURES)

Bernoulli function ``B(x)=\\frac{x}{e^x-1}`` for exponentially fitted upwinding.

The name `fbernoulli` has been chosen to avoid confusion with `Bernoulli` from JuliaStats/Distributions.jl

Returns a real number containing the result.

While `x/expm1(x)`  appears to be sufficiently accurate for all `x!=0`, 
it's derivative calculated via `ForwardDiff.jl` is not, 
so we use the polynomial approximation in the  
interval `(-bernoulli_small_threshold, bernoulli_small_threshold)`. 

Also, see the discussion in [#117](https://github.com/WIAS-PDELib/VoronoiFVM.jl/issues/117).
"""
function fbernoulli(x)
    return if x < -bernoulli_large_threshold
        -x
    elseif x > bernoulli_large_threshold
        zero(x)
    elseif abs(x) < bernoulli_small_threshold
        @inline bernoulli_horner(x)
    else
        x / expm1(x)
    end
end

##############################################################

""" 
$(SIGNATURES)

Bernoulli function ``B(x)=\\frac{x}{e^x-1}`` for exponentially
fitted upwind, joint evaluation for positive and negative
argument

Usually, we need ``B(x), B(-x)`` together, 
and it is cheaper to calculate them together.

Returns two real numbers containing the result for argument
`x` and argument `-x`. 

For the general approach to the implementation, see the discussion in [`fbernoulli`](@ref).
"""
function fbernoulli_pm(x)
    return if x < -bernoulli_large_threshold
        -x, zero(x)
    elseif x > bernoulli_large_threshold
        zero(x), x
    elseif abs(x) < bernoulli_small_threshold
        @inline y = bernoulli_horner(x)
        y, x + y
    else
        y = x / expm1(x)
        y, x + y
    end
end

"""
    $(SIGNATURES)
    
Non-pivoting inplace LU factorization using Doolittle's method.
Adapted from https://en.wikipedia.org/wiki/LU_decomposition#MATLAB_code_example.
"""
@inline function doolittle_ludecomp!(LU)
    n = size(LU, 1)
    @inbounds for i in 1:n
        for j in 1:(i - 1)
            for k in 1:(j - 1)
                LU[i, j] -= LU[i, k] * LU[k, j]
            end
            LU[i, j] /= LU[j, j]
        end
        for j in i:n
            for k in 1:(i - 1)
                LU[i, j] -= LU[i, k] * LU[k, j]
            end
        end
    end
    return LU
end

"""
$(SIGNATURES)

Non-pivoting inplace  upper and lower triangular solve of matrix
factorized with `doolittle_ludecomp!`.
Adapted from https://en.wikipedia.org/wiki/LU_decomposition#MATLAB_code_example.
"""
@inline function doolittle_lusolve!(LU, b)
    n = length(b)
    # LU= L+U-I
    # find solution of Ly = b and store it in b
    @inbounds for i in 1:n
        for k in 1:(i - 1)
            b[i] -= LU[i, k] * b[k]
        end
    end
    # find solution of Ux = b and store it in b
    @inbounds for i in n:-1:1
        for k in (i + 1):n
            b[i] -= LU[i, k] * b[k]
        end
        b[i] /= LU[i, i]
    end
    return b
end

"""
$(SIGNATURES)

Non-allocating, non-pivoting inplace solution of square linear system of equations `A*x=b`
using [Doolittle's method](https://en.wikipedia.org/wiki/LU_decomposition#Doolittle_algorithm).

After solution, `A` will contain the LU factorization, and `b` the result.

A pivoting version is available with Julia v1.9.
"""
@inline function inplace_linsolve!(A, b)
    doolittle_ludecomp!(A)
    return doolittle_lusolve!(A, b)
end

"""
$(SIGNATURES)

Non-allocating, pivoting, inplace solution of square linear system of equations `A*x=b`
using LU factorization from RecursiveFactorizations.jl.  

After solution, `A` will contain the LU factorization, and `b` the result.
`ipiv` must be an Int64 vector of the same length as `b`.
"""
@inline function inplace_linsolve!(A, b, ipiv)
    return LinearAlgebra.ldiv!(RecursiveFactorization.lu!(A, ipiv, Val(true), Val(false)), b)
end
