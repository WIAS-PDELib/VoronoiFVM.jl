#!/bin/sh
# -*-Julia-*-
#=
exec julia --startup-file=no --project= "$0" "$@"
=#

using Pkg
Pkg.activate(@__DIR__)

using Pluto

notebooks = [
    "nbproto.jl",
    "api-update.jl",
    "flux-reconstruction.jl",
    "problemcase.jl",
    "nonlinear-solvers.jl",
]


    # # if we add one package too many, this triggers the action of PlutoPkg
    # # when the NB is started
    # Pkg.add("AbstractTrees")
    # Pkg.update()

"""
    force_update_notebook_environment(notebook)

Forced update of notebook environment. `Pluto.update_notebook_environment(notebook)` heeds
the compat entries in the Project.toml and thus seems to do nothing.
"""
function force_update_notebook_environment(notebook)
    # cache the current environment
    thisproject=Pkg.project().path

    @info "Force updating packages in $(notebook):"

    # After this, notebook env + curent env are in sync
    Pluto.activate_notebook_environment(notebook)

    # Get list of current dependencies and their UUIDs:
    pkgs=[PackageSpec(uuid=v) for (k,v) in Pkg.project().dependencies]

    # Remove and re-add packages, thus ignoring compat
    Pkg.rm(pkgs)
    Pkg.add(pkgs)

    # let the environments sync
    sleep(1)

    # Sets compat to the current versions
    Pluto.update_notebook_environment(notebook)
    @info "Updating of  $(notebook) done\n"

    # Re-activate the current environment
    Pkg.activate(thisproject)
end

for notebook in notebooks
    force_update_notebook_environment(notebook)
end


dirs = ["test", "pluto-examples", "docs"]
for dir in dirs
    println("updating $(dir) environment")
    thisproject=Pkg.project().path
    Pkg.activate(joinpath(@__DIR__, "..", dir))
    Pkg.status()
    Pkg.update()
    Pkg.status()
    Pkg.activate(thisproject)
end
