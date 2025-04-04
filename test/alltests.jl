using ExplicitImports, Aqua
using ExampleJuggler: ExampleJuggler, cleanexamples, @testmodules, @testscripts
using VoronoiFVM: VoronoiFVM

ExampleJuggler.verbose!(true)
#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success. All files with filenames contained in
# the `ignore` list will not be included.
#
function run_tests_from_directory(testdir, prefix; ignore = [])
    @info "Directory $(testdir):"
    examples = filter(
        filename -> length(filename) >= length(prefix) &&
            filename[1:length(prefix)] == prefix &&
            filename ∉ ignore,
        basename.(readdir(testdir))
    )
    @info examples
    @testmodules(testdir, examples)
    return nothing
end

@testset "basictest" begin
    run_tests_from_directory(@__DIR__, "test_")
end

@testset "Development Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example0")
end
@testset "MultiD Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example5")
end
@testset "1D Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example1")
end
@testset "2D Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example2")
end

@testset "3D Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example3")
end

@testset "Misc Examples" begin
    run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example4")
end

# Run the notebooks as scripts in the test environment.
notebooks = [
    "nbproto.jl",
    "api-update.jl",
    "ode-diffusion1d.jl",
    "ode-wave1d.jl",
    "ode-nlstorage1d.jl",
    #"ode-brusselator.jl",
    "outflow.jl",
    "flux-reconstruction.jl",
    "interfaces1d.jl",
    "problemcase.jl",
    "nonlinear-solvers.jl",
    "heterogeneous-catalysis.jl",
]
@testset "Notebooks" begin
    @testscripts(joinpath(@__DIR__, "..", "pluto-examples"), notebooks)
end

@testset "ExplicitImports" begin
    @test ExplicitImports.check_no_implicit_imports(VoronoiFVM) === nothing
    @test ExplicitImports.check_no_stale_explicit_imports(VoronoiFVM) === nothing
end

@testset "Aqua" begin
    Aqua.test_ambiguities(VoronoiFVM, broken = true)
    Aqua.test_unbound_args(VoronoiFVM)
    Aqua.test_undefined_exports(VoronoiFVM)
    Aqua.test_project_extras(VoronoiFVM)
    Aqua.test_stale_deps(VoronoiFVM)
    Aqua.test_deps_compat(VoronoiFVM)
    Aqua.test_piracies(VoronoiFVM, broken = true)
    Aqua.test_persistent_tasks(VoronoiFVM)
end

@testset "UndocumentedNames" begin
    if isdefined(Docs, :undocumented_names) # >=1.11
        @test isempty(Docs.undocumented_names(VoronoiFVM))
    end
end
