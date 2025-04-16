# # 150: Impedance calculation
# ([source code](@__SOURCE_URL__))
#
#  Impedance calculation for
#
#    C u_t - (D u_x)_x + Ru = 0   in (0,1)
#      u(0,t)=1 + exp(iωt)
#      u(1,t)=0
#
#    Measurement: I(t)= D u_x(1,t)
#
#    Steady state:
#    - (D u0_x)_x + Ru0 = 0
#    u0(0,t)=1
#    u0(1,t)=0
#
#    Small signal ansatz for ω
#
#    u(x,t)= u0(x)+ ua(x) exp(iωt)
#
#    iωC ua - (D ua_x)_x + R u_a =0
#      ua(0)=1
#      ua(1)=0
#
#

module Example150_Impedance1D

using Printf
using VoronoiFVM
using ExtendableGrids: geomspace, simplexgrid
using GridVisualize

function main(;
        nref = 0, Plotter = nothing, verbose = false, unknown_storage = :sparse, assembly = :edgewise,
        L = 1.0, R = 1.0, D = 1.0, C = 1.0,
        ω0 = 1.0e-3, ω1 = 5.0e1
    )

    # Create array which is refined close to 0
    h0 = 0.005 / 2.0^nref
    h1 = 0.1 / 2.0^nref

    X = geomspace(0, L, h0, h1)

    # Create discretitzation grid
    grid = simplexgrid(X)

    # Create and fill data
    data = (R = R, D = D, C = C)

    # Declare constitutive functions
    flux = function (f, u, edge, data)
        f[1] = data.D * (u[1, 1] - u[1, 2])
        return nothing
    end

    storage = function (f, u, node, data)
        f[1] = data.C * u[1]
        return nothing
    end

    reaction = function (f, u, node, data)
        f[1] = data.R * u[1]
        return nothing
    end

    excited_bc = 1
    excited_bcval = 1.0
    excited_spec = 1
    meas_bc = 2

    # Create physics struct
    physics = VoronoiFVM.Physics(;
        data = data,
        flux = flux,
        storage = storage,
        reaction = reaction
    )
    # Create discrete system and enable species
    sys = VoronoiFVM.System(grid, physics; unknown_storage = unknown_storage, assembly = assembly)

    enable_species!(sys, excited_spec, [1])

    # Create test functions for current measurement

    factory = TestFunctionFactory(sys)
    measurement_testfunction = testfunction(factory, [excited_bc], [meas_bc])

    boundary_dirichlet!(sys, excited_spec, excited_bc, excited_bcval)
    boundary_dirichlet!(sys, excited_spec, meas_bc, 0.0)

    steadystate = solve(sys)

    function meas_stdy(meas, U)
        u = reshape(U, sys)
        meas[1] = -VoronoiFVM.integrate_stdy(sys, measurement_testfunction, u)[excited_spec]
        return nothing
    end

    function meas_tran(meas, U)
        u = reshape(U, sys)
        meas[1] = -VoronoiFVM.integrate_tran(sys, measurement_testfunction, u)[excited_spec]
        return nothing
    end

    dmeas_stdy = measurement_derivative(sys, meas_stdy, steadystate)
    dmeas_tran = measurement_derivative(sys, meas_tran, steadystate)

    # Create Impeadancs system from steady state
    isys = VoronoiFVM.ImpedanceSystem(sys, steadystate, excited_spec, excited_bc)

    # Prepare recording of impedance results
    allomega = zeros(0)

    # for calculated data
    allI0 = zeros(Complex{Float64}, 0)
    allIL = zeros(Complex{Float64}, 0)

    # for exact data
    allIx0 = zeros(Complex{Float64}, 0)
    allIxL = zeros(Complex{Float64}, 0)

    ω = ω0

    UZ = unknowns(isys)
    while ω < ω1

        # solve impedance system
        solve!(UZ, isys, ω)

        # calculate approximate solution
        # obtain measurement in frequency  domain
        IL = impedance(isys, ω, steadystate, dmeas_stdy, dmeas_tran)

        # record approximate solution
        push!(allomega, ω)
        push!(allIL, IL)

        # record exact solution
        iω = 1im * ω
        z = sqrt(iω * data.C / data.D + data.R / data.D)
        eplus = exp(z * L)
        eminus = exp(-z * L)
        IxL = 2.0 * data.D * z / (eplus - eminus)

        push!(allIxL, 1 / IxL)

        # increase omega
        ω = ω * 1.1
    end

    p = GridVisualizer(; Plotter = Plotter)
    scalarplot!(
        p, real(allIxL), imag(allIxL); label = "exact", color = :red,
        linestyle = :dot
    )
    scalarplot!(
        p, real(allIL), imag(allIL); label = "calc", show = true, clear = false,
        color = :blue, linestyle = :solid
    )

    return sum(allIL)
end

using Test
function runtests()
    testval = 57.92710286186797 + 23.163945443946027im
    @test main(; unknown_storage = :sparse, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :edgewise) ≈ testval &&
        main(; unknown_storage = :sparse, assembly = :cellwise) ≈ testval &&
        main(; unknown_storage = :dense, assembly = :cellwise) ≈ testval
    return nothing
end

end
