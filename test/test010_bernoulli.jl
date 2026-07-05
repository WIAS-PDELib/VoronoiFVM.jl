module test_bernoulli
using VoronoiFVM
using Test

function runtests()
    function B_Big(x)
        bx = BigFloat(x)
        return Float64(bx / (exp(bx) - one(bx)))
    end
    smallrange = -1:1.00001e-5:1
    largerange = -100:1.00001e-3:100
    maxerror(X, b, f) = maximum(abs.(b.(X) .- f.(X)))
    @test maxerror(smallrange, B_Big, fbernoulli) < 1.0e-14
    @test maxerror(largerange, B_Big, fbernoulli) < 1.0e-14

    @test maxerror(smallrange, B_Big, (x) -> fbernoulli_pm(x)[1]) < 1.0e-14
    @test maxerror(largerange, B_Big, (x) -> fbernoulli_pm(x)[1]) < 1.0e-14

    @test maxerror(smallrange, x -> B_Big(-x), (x) -> fbernoulli_pm(x)[2]) < 1.0e-14
    @test maxerror(largerange, x -> B_Big(-x), (x) -> fbernoulli_pm(x)[2]) < 1.0e-14

    return true
end

end
