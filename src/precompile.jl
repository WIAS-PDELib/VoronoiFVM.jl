## Attempt on https://discourse.julialang.org/t/22-seconds-to-3-and-now-more-lets-fix-all-of-the-differentialequations-jl-universe-compile-times/66313
using PrecompileTools
@static if VERSION > v"1.7" && !Sys.isapple()
    PrecompileTools.@compile_workload let
        function lin1()
            n = 5
            X = 0:(1.0 / n):1
            grid = simplexgrid(X, X)
            ispec = 1

            function diffusion_flux(D::T) where {T}
                (y, u, edge) -> y[1] = D(u[1, 1] + u[1, 2]) * (u[1, 1] - u[1, 2])
            end

            flux!(y, u, edge) = y[1] = u[1, 1] - u[1, 2]

            reaction!(y, u, node) = y[1] = u[1, 1]^3 - 10

            function bc!(args...)
                boundary_dirichlet!(args...; region = 1, value = 0)
                boundary_dirichlet!(args...; region = 3, value = 1)
            end

            solve(
                VoronoiFVM.System(
                    simplexgrid(X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :cellwise
                ); verbose = ""
            )
            solve(
                VoronoiFVM.System(
                    simplexgrid(X, X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :cellwise
                ); verbose = ""
            )
            solve(
                VoronoiFVM.System(
                    simplexgrid(X, X, X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :cellwise
                ); verbose = ""
            )
            solve(
                VoronoiFVM.System(
                    simplexgrid(X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :edgewise
                ); verbose = ""
            )
            solve(
                VoronoiFVM.System(
                    simplexgrid(X, X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :edgewise
                ); verbose = ""
            )
            solve(
                VoronoiFVM.System(
                    simplexgrid(X, X, X); species = 1, flux = flux!, reaction = reaction!, breaction = bc!,
                    assembly = :edgewise
                ); verbose = ""
            )
        end
        lin1()
    end
end
