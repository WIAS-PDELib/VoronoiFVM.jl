# Development hints
Here, some development hints are given which mainly concern tests and documentation generation.

## ReTest
[ReTest.jl](https://github.com/JuliaTesting/ReTest.jl) allows to run tests from a subset of testsets specified by regular expression. For the gemeral workflow, see the [quick start](https://juliatesting.github.io/ReTest.jl/stable/#Quick-start) documentation.

Testing the package code can be done just by calling `test/runtests.jl`. This file is also invoked by the CI, and it just includes `test/alltests.jl`. 

Alternatively, test code can be run via ReTest. For this purpose, the module `test/VoronoiFVMReTest.jl` can be loaded in an environment with `VoronoiFVM.jl` (developed), test dependencies (see `test/Project.toml`), `ReTest.jl` and `Revise.jl`. A shared environment (e.g. via `julia --project=@VoronoiFVM`) is a convenient way to maintain such an environment locally.

After loading `VoronoiFVM` via `ReTest` as follows
```
julia> using ReTest
julia> ReTest.load("test/VoronoiFVMReTest.jl")
```
e.g.
```
julia> retest("Aqua")
```
runs just the "Aqua" testset. Via
```
julia> retest(dry=true)
```
one obtains a list of all possible testsets.


## Pluto notebooks
The pluto notebooks in this package are "triple use":
- As typical Pluto notebooks, they are self-contained in the senses that they contain their own Project and Manifest files. So users can just download and execute them.
- If they run with the environment variable `PLUTO_PROJECT` set to some julia environment, this environment will activated at the start of the notebook. In particular, they use `Revise.jl` so they can be run during development of VoronoiFVM.jl. See also  https://github.com/fonsp/Pluto.jl/issues/1788 .
- During CI tests (including ReTest), they are run as scripts. For this purpose they are wrapped into temporary modules, and @test macros can be used in the notebooks.
