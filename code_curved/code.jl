# Install dependencies
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Load dependencies
using LinearAlgebra
using Statistics
using DelimitedFiles: DelimitedFiles, writedlm, readdlm

using Trixi
using OrdinaryDiffEq
using DiffEqCallbacks
using SummationByPartsOperators


import PolyesterWeave, ThreadingUtilities

using LaTeXStrings
using Plots: Plots, plot, plot!, scatter, scatter!, savefig

using Trixi2Vtk: trixi2vtk

using PrettyTables: PrettyTables, pretty_table, ft_printf


const figdir = joinpath(dirname(@__DIR__), "figures")


function plot_kwargs()
    fontsizes = (
        xtickfontsize = 14, ytickfontsize = 14,
        xguidefontsize = 16, yguidefontsize = 16,
        legendfontsize = 14)
    (; linewidth = 3, gridlinewidth = 2,
        markersize = 8, markerstrokewidth = 4,
        fontsizes..., size=(600, 500))
end

################################################################################
# Convergence tests

function convergence_tests_curvilinear_2d_euler(; latex = false)
    # Directly test the 4-2, 6-3, and 8-4 upwind operators; EOCs for others are similar
    orders = [4, 6, 8]

    for accuracy_order in orders
        if accuracy_order == 8
            mesh_file_name = "mesh_convergence_test_deg2.mesh"
        else
            mesh_file_name = "mesh_convergence_test_deg1.mesh"
        end

        num_nodes = [17, 34, 68, 136, 272]
        splitting = splitting_vanleer_haenel
        @info "2D compressible Euler equations" accuracy_order splitting
        _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
                                                accuracy_order, splitting, latex)

        num_nodes = [17, 34, 68, 136, 272]
        splitting = splitting_drikakis_tsangaris
        @info "2D compressible Euler equations" accuracy_order splitting
        _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
                                                accuracy_order, splitting, latex)

        num_nodes = [17, 34, 68, 136, 272]
        splitting = splitting_lax_friedrichs
        @info "2D compressible Euler equations" accuracy_order splitting
        _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
                                                accuracy_order, splitting, latex)
    end

    return nothing
end

function _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes, accuracy_order,
                                                   splitting, latex = false)
    num_elements = Vector{Int}()
    errors = Vector{Float64}()

    for nnodes in num_nodes
        nelements = 16 # both meshes have fixed 16 elements so this value is just hard coded
        tol = 1.0e-13
        res = compute_errors_curvilinear_2d_euler(; mesh_file_name, nnodes,
                                                  accuracy_order, tol, splitting)
        push!(num_elements, nelements)
        push!(errors, first(res.l2))
    end

    eoc = compute_eoc(num_nodes, errors)

    # print results;
    data = hcat(num_elements, num_nodes, errors, eoc)
    header = ["#elements", "#nodes", "L2 error", "L2 EOC"]
    kwargs = (; header, formatters=(ft_printf("%3d", [1, 2]),
                                    ft_printf("%.2e", [3]),
                                    ft_printf("%.2f", [4])))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end

    return nothing
end

function compute_errors_curvilinear_2d_euler(; mesh_file_name, nnodes,
                                             accuracy_order, tol, splitting)
    equations = CompressibleEulerEquations2D(1.4)

    # Modify the manufactured solution test to use `L = sqrt(2)` in the initial and source terms
    function initial_condition_convergence_shifted(x, t, equations::CompressibleEulerEquations2D)
        c = 2
        A = 0.1
        L = sqrt(2)
        f = 1 / L
        ω = 2 * pi * f
        ini = c + A * sin(ω * (x[1] + x[2] - t))

        rho = ini
        rho_v1 = ini
        rho_v2 = ini
        rho_e = ini^2

        return SVector(rho, rho_v1, rho_v2, rho_e)
    end

    @inline function source_terms_convergence_shifted(u, x, t,
                                                      equations::CompressibleEulerEquations2D)
        # Same settings as in `initial_condition`
        c = 2
        A = 0.1
        L = sqrt(2)
        f = 1 / L
        ω = 2 * pi * f
        γ = equations.gamma

        x1, x2 = x
        si, co = sincos(ω * (x1 + x2 - t))
        rho = c + A * si
        rho_x = ω * A * co
        # Note that d/dt rho = -d/dx rho = -d/dy rho.

        tmp = (2 * rho - 1) * (γ - 1)

        du1 = rho_x
        du2 = rho_x * (1 + tmp)
        du3 = du2
        du4 = 2 * rho_x * (rho + tmp)

        return SVector(du1, du2, du3, du4)
    end
    initial_condition = initial_condition_convergence_shifted
    source_terms = source_terms_convergence_shifted

    D_upw = upwind_operators(
        SummationByPartsOperators.Mattsson2017;
        derivative_order = 1,
        accuracy_order,
        xmin = -1.0, xmax = 1.0,
        N = nnodes)
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralStrongForm(FluxUpwind(splitting)),
                   volume_integral = VolumeIntegralUpwind(splitting))

    mesh_file = joinpath(@__DIR__, mesh_file_name)

    mesh = UnstructuredMesh2D(mesh_file, periodicity = true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
                                        solver; source_terms)

    ode = semidiscretize(semi, (0.0, 2.0))
    sol = solve(ode, SSPRK43(); ode_default_options()...,
                abstol = tol, reltol = tol)

    analysis_callback = AnalysisCallback(semi)
    return analysis_callback(sol)
end

function compute_eoc(Ns, errors)
    eoc = similar(errors)
    eoc[begin] = NaN # no EOC defined for the first grid
    for idx in Iterators.drop(eachindex(errors, Ns, eoc), 1)
        eoc[idx] = -( log(errors[idx] / errors[idx - 1]) / log(Ns[idx] / Ns[idx - 1]) )
    end
    return eoc
end

################################################################################
# Free-stream preservation tests

function free_stream_preservation_generate_data()
    for flux_splitting in [splitting_lax_friedrichs, splitting_vanleer_haenel, splitting_drikakis_tsangaris]
        for accuracy_order in 2:9
            @info "Generate FSP data for linear boundaries with" accuracy_order flux_splitting
            t, error_density = compute_fsp_error_degOne(; accuracy_order, flux_splitting)
            open(joinpath(figdir, "fsp_error_degOne_$(accuracy_order)_"*string(flux_splitting)*".dat"), "w") do io
                println(io, "# t\tL2_error_density")
                writedlm(io, hcat(t, error_density))
            end

            @info "Generate FSP data for quadratic boundaries with" accuracy_order flux_splitting
            t, error_density = compute_fsp_error_degTwo(; accuracy_order, flux_splitting)
            open(joinpath(figdir, "fsp_error_degTwo_$(accuracy_order)_"*string(flux_splitting)*".dat"), "w") do io
                println(io, "# t\tL2_error_density")
                writedlm(io, hcat(t, error_density))
            end
        end
    end
end

# Note the default tolerances for the FSP testing are taken to be quite strict
function compute_fsp_error_degOne(; accuracy_order = 4, nnodes = 17,
                                  flux_splitting = splitting_vanleer_haenel,
                                  source_of_coefficients = Mattsson2017,
                                  tspan = (0.0, 10.0),
                                  tol = 1.0e-12)
    equations = CompressibleEulerEquations2D(1.4)

    initial_condition = initial_condition_constant

    # Boundary conditions for free-stream preservation test
    boundary_condition_free_stream = BoundaryConditionDirichlet(initial_condition)

    boundary_conditions = Dict(:outerCircle => boundary_condition_free_stream,
                               :cone1 => boundary_condition_free_stream,
                               :cone2 => boundary_condition_free_stream,
                               :iceCream => boundary_condition_free_stream)

    # Use upwind SBP discretization
    D_upw = upwind_operators(source_of_coefficients;
                             derivative_order = 1,
                             accuracy_order,
                             xmin = -1.0, xmax = 1.0,
                             N = nnodes)
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralStrongForm(FluxUpwind(flux_splitting)),
                   volume_integral = VolumeIntegralUpwind(flux_splitting))

    # unstructured mesh with bi-linear elements
    mesh_file = joinpath(@__DIR__, "mesh_inner_outer_boundaries_deg1.mesh")

    mesh = UnstructuredMesh2D(mesh_file)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                        boundary_conditions=boundary_conditions)

    ode = semidiscretize(semi, tspan)

    saveat = tspan[end]
    saved_values = SavedValues(Float64, Float64)
    save_func = let cb = AnalysisCallback(semi)
        function save_func(u_ode, t, integrator)
            semi = integrator.p
            analysis_callback = cb.affect!
            (; analyzer) = analysis_callback
            cache_analysis = analysis_callback.cache

            l2_error, linf_error = Trixi.calc_error_norms(u_ode, t,
                                                          analyzer, semi,
                                                          cache_analysis)
            return first(l2_error)
        end
    end
    saving = SavingCallback(save_func, saved_values; saveat)

    sol = solve(ode, SSPRK43();
                abstol = tol, reltol = tol,
                ode_default_options()..., callback = saving, tstops = saveat)
    return (; t = saved_values.t, error_density = saved_values.saveval)
end

# Note the default tolerances for the FSP testing are taken to be quite strict
function compute_fsp_error_degTwo(; accuracy_order = 4, nnodes = 17,
                                  flux_splitting = splitting_vanleer_haenel,
                                  source_of_coefficients = Mattsson2017,
                                  tspan = (0.0, 10.0),
                                  tol = 1.0e-12)
    equations = CompressibleEulerEquations2D(1.4)

    initial_condition = initial_condition_constant

    # Boundary conditions for free-stream preservation test
    boundary_condition_free_stream = BoundaryConditionDirichlet(initial_condition)

    boundary_conditions = Dict(:outerCircle => boundary_condition_free_stream,
                               :cone1 => boundary_condition_free_stream,
                               :cone2 => boundary_condition_free_stream,
                               :iceCream => boundary_condition_free_stream)

    # Use upwind SBP discretization
    D_upw = upwind_operators(source_of_coefficients;
                             derivative_order = 1,
                             accuracy_order,
                             xmin = -1.0, xmax = 1.0,
                             N = nnodes)
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralStrongForm(FluxUpwind(flux_splitting)),
                   volume_integral = VolumeIntegralUpwind(flux_splitting))

    # unstructured mesh with quadratic elements
    mesh_file = joinpath(@__DIR__, "mesh_inner_outer_boundaries_deg2.mesh")

    mesh = UnstructuredMesh2D(mesh_file)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                        boundary_conditions=boundary_conditions)

    ode = semidiscretize(semi, tspan)

    saveat = tspan[end]
    saved_values = SavedValues(Float64, Float64)
    save_func = let cb = AnalysisCallback(semi)
        function save_func(u_ode, t, integrator)
            semi = integrator.p
            analysis_callback = cb.affect!
            (; analyzer) = analysis_callback
            cache_analysis = analysis_callback.cache

            l2_error, linf_error = Trixi.calc_error_norms(u_ode, t,
                                                          analyzer, semi,
                                                          cache_analysis)
            return first(l2_error)
        end
    end
    saving = SavingCallback(save_func, saved_values; saveat)

    sol = solve(ode, SSPRK43();
                abstol = tol, reltol = tol,
                ode_default_options()..., callback = saving, tstops = saveat)
    return (; t = saved_values.t, error_density = saved_values.saveval)
end

################################################################################
# Isentropic vortex on warped mesh

function isentropic_vortex_generate_data()
    for flux_splitting in [splitting_lax_friedrichs, splitting_drikakis_tsangaris, splitting_vanleer_haenel]
        for accuracy_order in 8:8 # only use the 8-4 operator for this curved mesh
            @info "Generate data for" accuracy_order flux_splitting
            t, error_density = compute_error_isentropic_vortex(; accuracy_order,
                                                                 nnodes = 17,
                                                                 source_of_coefficients = Mattsson2017
                                                                 )
            open(joinpath(figdir, "isentropic_vortex_order_$(accuracy_order)_"*string(flux_splitting)*"_warped.dat"), "w") do io
                println(io, "# t\tL2_error_density")
                writedlm(io, hcat(t, error_density))
            end
        end
    end
end

function isentropic_vortex_plot_results()
    fig = plot(xguide = L"Time $t$", yguide = L"$L^2$ error of the density";
               xscale = :log10, yscale = :log10,
               plot_kwargs()...)

    linestyles = [:solid, :dash, :dot]
    for(splitting_style, linestyle) in zip([splitting_lax_friedrichs, splitting_drikakis_tsangaris, splitting_vanleer_haenel], linestyles)
        data = readdlm(joinpath(figdir,  "isentropic_vortex_order_8_"*string(splitting_style)*"_warped.dat"), comments = true)
        plot!(fig, data[:, 1], data[:, 2];
              label = splitting_style, linestyle,
              plot_kwargs()...)
    end

    plot!(fig, legend = :bottomright)
    savefig(fig, joinpath(figdir, "isentropic_vortex_error_warped.pdf"))
    @info "Error plot saved in the directory `figdir`" figdir
    return nothing
end

function compute_error_isentropic_vortex(; accuracy_order = 8, nnodes = 17,
                                           flux_splitting = splitting_vanleer_haenel,
                                           source_of_coefficients = Mattsson2017,
                                           tspan = (0.0, 1000.0),
                                           tol = 1.0e-6)
    equations = CompressibleEulerEquations2D(1.4)

    """
        initial_condition_isentropic_vortex(x, t, equations)

    The classical isentropic vortex test case of
    - Chi-Wang Shu (1997)
    Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory
    Schemes for Hyperbolic Conservation Laws.
    [NASA/CR-97-206253](https://ntrs.nasa.gov/citations/19980007543)
    """
    function initial_condition(x, t, equations::CompressibleEulerEquations2D)
        # needs appropriate mesh size, e.g. [-10,-10]x[10,10]
        # for error convergence: make sure that the end time is such that the is back at the initial state!!
        # for the current velocity and domain size: t_end should be a multiple of 20s
        # initial center of the vortex
        inicenter = SVector(0.0, 0.0)
        # size and strength of the vortex
        iniamplitude = 5.0
        # base flow
        rho = 1.0
        v1 = 1.0
        v2 = 1.0
        vel = SVector(v1, v2)
        p = 25.0
        rt = p / rho                  # ideal gas equation
        t_loc = 0.0
        cent = inicenter + vel * t_loc      # advection of center
        # ATTENTION: handle periodic BC, but only for v1 = v2 = 1.0 (!!!!)

        cent = x - cent # distance to center point

        #cent=cross(iniaxis,cent) # distance to axis, tangent vector, length r
        # cross product with iniaxis = [0, 0, 1]
        cent = SVector(-cent[2], cent[1])
        r2 = cent[1]^2 + cent[2]^2
        du = iniamplitude / (2 * π) * exp(0.5 * (1 - r2)) # vel. perturbation
        dtemp = -(equations.gamma - 1) / (2 * equations.gamma * rt) * du^2 # isentropic
        rho = rho * (1 + dtemp)^(1 / (equations.gamma - 1))
        vel = vel + du * cent
        v1, v2 = vel
        p = p * (1 + dtemp)^(equations.gamma / (equations.gamma - 1))
        prim = SVector(rho, v1, v2, p)
        return prim2cons(prim, equations)
    end

    # Use upwind SBP discretization
    D_upw = upwind_operators(source_of_coefficients;
                            derivative_order = 1,
                            accuracy_order,
                            xmin = -1.0, xmax = 1.0,
                            N = nnodes)
    solver = FDSBP(D_upw,
                    surface_integral = SurfaceIntegralStrongForm(FluxUpwind(flux_splitting)),
                    volume_integral = VolumeIntegralUpwind(flux_splitting))

    mesh_file = joinpath(@__DIR__, "mesh_warped_vortex.mesh")

    mesh = UnstructuredMesh2D(mesh_file, periodicity=true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    ode = semidiscretize(semi, tspan)

    saveat = range(tspan..., step = 20)[2:end]
    saved_values = SavedValues(Float64, Float64)
    save_func = let cb = AnalysisCallback(semi)
        function save_func(u_ode, t, integrator)
            semi = integrator.p
            analysis_callback = cb.affect!
            (; analyzer) = analysis_callback
            cache_analysis = analysis_callback.cache

            l2_error, linf_error = Trixi.calc_error_norms(u_ode, t,
                                                          analyzer, semi,
                                                          cache_analysis)
            return first(l2_error)
        end
    end
    saving = SavingCallback(save_func, saved_values; saveat)

    sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True());
                abstol = tol, reltol = tol, maxiters = 1e10,
                ode_default_options()..., callback = saving, tstops = saveat)
    return (; t = saved_values.t, error_density = saved_values.saveval)
end

# Used to create figure data for the ParaView plots where the h5 files
# are saved to the folder "out"
function run_isentropic_vortex(; accuracy_order = 8, nnodes = 17,
                                 flux_splitting = splitting_vanleer_haenel,
                                 source_of_coefficients = Mattsson2017,
                                 tspan = (0.0, 20.0),
                                 tol = 1.0e-6)
    equations = CompressibleEulerEquations2D(1.4)

    """
        initial_condition_isentropic_vortex(x, t, equations)

    The classical isentropic vortex test case of
    - Chi-Wang Shu (1997)
    Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory
    Schemes for Hyperbolic Conservation Laws.
    [NASA/CR-97-206253](https://ntrs.nasa.gov/citations/19980007543)
    """
    function initial_condition(x, t, equations::CompressibleEulerEquations2D)
        # needs appropriate mesh size, e.g. [-10,-10]x[10,10]
        # for error convergence: make sure that the end time is such that the is back at the initial state!!
        # for the current velocity and domain size: t_end should be a multiple of 20s
        # initial center of the vortex
        inicenter = SVector(0.0, 0.0)
        # size and strength of the vortex
        iniamplitude = 5.0
        # base flow
        rho = 1.0
        v1 = 1.0
        v2 = 1.0
        vel = SVector(v1, v2)
        p = 25.0
        rt = p / rho                  # ideal gas equation
        t_loc = 0.0
        cent = inicenter + vel * t_loc      # advection of center
        # ATTENTION: handle periodic BC, but only for v1 = v2 = 1.0 (!!!!)

        cent = x - cent # distance to center point

        #cent=cross(iniaxis,cent) # distance to axis, tangent vector, length r
        # cross product with iniaxis = [0, 0, 1]
        cent = SVector(-cent[2], cent[1])
        r2 = cent[1]^2 + cent[2]^2
        du = iniamplitude / (2 * π) * exp(0.5 * (1 - r2)) # vel. perturbation
        dtemp = -(equations.gamma - 1) / (2 * equations.gamma * rt) * du^2 # isentropic
        rho = rho * (1 + dtemp)^(1 / (equations.gamma - 1))
        vel = vel + du * cent
        v1, v2 = vel
        p = p * (1 + dtemp)^(equations.gamma / (equations.gamma - 1))
        prim = SVector(rho, v1, v2, p)
        return prim2cons(prim, equations)
    end

    # Use upwind SBP discretization
    D_upw = upwind_operators(source_of_coefficients;
                            derivative_order = 1,
                            accuracy_order,
                            xmin = -1.0, xmax = 1.0,
                            N = nnodes)
    solver = FDSBP(D_upw,
                    surface_integral = SurfaceIntegralStrongForm(FluxUpwind(flux_splitting)),
                    volume_integral = VolumeIntegralUpwind(flux_splitting))

    mesh_file = joinpath(@__DIR__, "mesh_warped_vortex.mesh")

    mesh = UnstructuredMesh2D(mesh_file, periodicity=true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 1000
    analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

    alive_callback = AliveCallback(analysis_interval = analysis_interval)

    save_solution = SaveSolutionCallback(dt = 5.0,
                                         save_initial_solution = true,
                                         save_final_solution = true,
                                         solution_variables = cons2prim)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback,
                            save_solution,
                            alive_callback)

    sol = solve(ode, SSPRK43(); abstol = 1.0e-6, reltol = 1.0e-6, dt = 1e-3,
                ode_default_options()..., callback = callbacks)
    summary_callback()
end