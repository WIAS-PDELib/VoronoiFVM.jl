using Documenter, ExampleJuggler, PlutoStaticHTML, VoronoiFVM, DocumenterCitations
using ExtendableGrids, GridVisualize, LinearAlgebra, RecursiveArrayTools, SciMLBase

using ExtendableFEMBase, ExtendableFEM

using OrdinaryDiffEqBDF, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqTsit5

function make(;
        with_examples = true,
        with_notebooks = true
    )

    bib = CitationBibliography(
        joinpath(@__DIR__, "..", "CITATIONS.bib");
        style = :numeric
    )

    ExampleJuggler.verbose!(true)

    cleanexamples()
    notebookdir = joinpath(@__DIR__, "..", "pluto-examples")
    exampledir = joinpath(@__DIR__, "..", "examples")

    size_threshold_ignore = []

    pages = [
        "Home" => "index.md",
        "Changelog" => "changes.md",
        "method.md",
        "API Documentation" => [
            "system.md",
            "physics.md",
            "solutions.md",
            "solver.md",
            "post.md",
            "quantities.md",
            "misc.md",
            "internal.md",
            "allindex.md",
            "devel.md",
            "extensions.md",
        ],
    ]


    if with_notebooks
        notebooks = [
            "OrdinaryDiffEq.jl nonlinear diffusion" => "ode-diffusion1d.jl",
            "OrdinaryDiffEq.jl 1D wave equation" => "ode-wave1d.jl",
            "OrdinaryDiffEq.jl changing mass matrix" => "ode-nlstorage1d.jl",
            "OrdinaryDiffEq.jl brusselator" => "ode-brusselator.jl",
            "Coupling with Catalyst.jl" => "heterogeneous-catalysis.jl",
            "Outflow boundary conditions" => "outflow.jl",
            "Obtaining vector fields" => "flux-reconstruction.jl",
            "Internal interfaces (1D)" => "interfaces1d.jl",
            "A case for caution" => "problemcase.jl",
            "Nonlinear solver control" => "nonlinear-solvers.jl",
            "Bernoulli function test" => "bernoulli.jl",
            "API Updates" => "api-update.jl",
        ]
        notebook_examples = @docplutonotebooks(notebookdir, notebooks, iframe = false)
        notebook_examples = vcat(["About the notebooks" => "notebooks.md"], notebook_examples)
        size_threshold_ignore = last.(notebook_examples)
        push!(pages, "Tutorial Notebooks" => notebook_examples)
    end

    if with_examples
        modules = filter(ex -> splitext(ex)[2] == ".jl", basename.(readdir(exampledir)))
        module_examples = @docmodules(exampledir, modules, use_module_titles = true)
        module_examples = vcat(["About the examples" => "runexamples.md"], module_examples)
        push!(pages, "Examples" => module_examples)
    end

    makedocs(;
        sitename = "VoronoiFVM.jl",
        modules = [
            VoronoiFVM,
            # define extension modules manually: https://github.com/JuliaDocs/Documenter.jl/issues/2124#issuecomment-1557473415
            Base.get_extension(VoronoiFVM, :VoronoiFVMExtendableFEMBaseExt),
        ],
        plugins = [bib],
        checkdocs = :all,
        clean = false,
        doctest = false,
        authors = "J. Fuhrmann",
        repo = "https://github.com/WIAS-PDELib/VoronoiFVM.jl",
        format = Documenter.HTML(;
            size_threshold_ignore,
            assets = String["assets/citations.css"],
            mathengine = MathJax3()
        ),
        pages
    )


    cleanexamples()

    return if !isinteractive()
        deploydocs(; repo = "github.com/WIAS-PDELib/VoronoiFVM.jl.git")
    end
end

if isinteractive()
    make(; with_examples = false, with_notebooks = false)
else
    make(; with_examples = true, with_notebooks = true)
end
