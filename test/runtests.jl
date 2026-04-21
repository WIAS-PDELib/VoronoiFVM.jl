using Test
#include("alltests.jl")

using ExampleJuggler: ExampleJuggler, @testmodules
ExampleJuggler.verbose!(true)
@testmodules(joinpath(@__DIR__, "..", "examples"), ["Example440_ParallelState.jl"])
