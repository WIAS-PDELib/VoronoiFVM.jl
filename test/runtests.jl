using Test, ExplicitImports, Aqua
using ExampleJuggler: ExampleJuggler, cleanexamples, @testmodules, @testscripts
using VoronoiFVM: VoronoiFVM

ExampleJuggler.verbose!(true)
#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success.
#
function run_tests_from_directory(testdir, prefix)
    @info "Directory $(testdir):"
    examples = filter(ex -> length(ex) >= length(prefix) && ex[1:length(prefix)] == prefix, basename.(readdir(testdir)))
    @info examples
    @testmodules(testdir, examples)
    return nothing
end

function run_all_tests(; run_notebooks = false, notebooksonly = false)
    if !notebooksonly
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
    end

    if run_notebooks
        ENV["VORONOIFVM_RUNTESTS"] = "1"
        notebooks = [
            "nbproto.jl",
            "ode-diffusion1d.jl",
            "ode-wave1d.jl",
            "ode-nlstorage1d.jl",
            "ode-brusselator.jl",
            "outflow.jl",
            "flux-reconstruction.jl",
            "interfaces1d.jl",
            "problemcase.jl",
            "nonlinear-solvers.jl",
            "api-update.jl",
            "heterogeneous-catalysis.jl",
        ]
        @testset "Notebooks" begin
            @testscripts(joinpath(@__DIR__, "..", "pluto-examples"), notebooks)
        end
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

    if isdefined(Docs, :undocumented_names) # >=1.11
        @testset "UndocumentedNames" begin
            @test isempty(Docs.undocumented_names(VoronoiFVM))
        end
    end
    return nothing
end

# Don't run notebooks on 1.12: https://github.com/fonsp/Pluto.jl/issues/2939
if haskey(ENV, "EXAMPLES_ONLY")
    run_all_tests(; run_notebooks = false, notebooksonly = false)
elseif haskey(ENV, "NOTEBOOKS_ONLY")
    run_all_tests(; run_notebooks = true, notebooksonly = true)
else
    run_all_tests(; run_notebooks = VERSION < v"1.12.0-DEV.0", notebooksonly = false)
end
