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


################################################################################
# Spectra

function plot_spectra()
    # Common axis formatting
    xlims = (-320, 0)
    ylims = (-250, 250)

    # Show spectra for comparable upwind SBP FD and DGSEM schemes in a single plot
    let accuracy_order = 4, polydeg = 2
        fig = plot(xguide = "Real part", yguide = "Imaginary part", xlims=xlims, ylims=ylims)
        let
            nnodes = 48
            initial_refinement_level = 2
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :solid)
        end
        let
            nnodes = 24
            initial_refinement_level = 3
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dash)
        end
        let
            nnodes = 12
            initial_refinement_level = 4
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dot)
        end
        let
            initial_refinement_level = 6
            λ = compute_spectrum_1d_dgsem(; initial_refinement_level, polydeg)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) DG elements with polynomial degree $polydeg"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dashdotdot,
                  color = :black)
        end
        plot!(fig, legend = :outertop)
        savefig(fig, joinpath(figdir,
                "spectra_linear_advection_1d_order$(accuracy_order)_polydeg$(polydeg).pdf"))
    end

    let accuracy_order = 6, polydeg = 3
        fig = plot(xguide = "Real part", yguide = "Imaginary part", xlims=xlims, ylims=ylims)
        let
            nnodes = 64
            initial_refinement_level = 2
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :solid)
        end
        let
            nnodes = 32
            initial_refinement_level = 3
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dash)
        end
        let
            nnodes = 16
            initial_refinement_level = 4
            λ = compute_spectrum_1d(; initial_refinement_level, nnodes, accuracy_order)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) upwind SBP elements with $(nnodes) nodes"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dot)
        end
        let
            initial_refinement_level = 6
            λ = compute_spectrum_1d_dgsem(; initial_refinement_level, polydeg)
            @show extrema(real, λ)
            @show maximum(abs, λ)
            λ = sort_spectrum(λ)
            label = "$(2^initial_refinement_level) DG elements with polynomial degree $polydeg"
            plot!(fig, real.(λ), imag.(λ); label, plot_kwargs()..., linestyle = :dashdotdot,
                  color = :black)
        end
        plot!(fig, legend = :outertop)
        savefig(fig, joinpath(figdir,
                "spectra_linear_advection_1d_order$(accuracy_order)_polydeg$(polydeg).pdf"))
    end

    @info "1D spectra saved in the directory `figdir`" figdir
    return nothing
end

function plot_kwargs()
    fontsizes = (
        xtickfontsize = 14, ytickfontsize = 14,
        xguidefontsize = 16, yguidefontsize = 16,
        legendfontsize = 14)
    (; linewidth = 3, gridlinewidth = 2,
        markersize = 8, markerstrokewidth = 4,
        fontsizes..., size=(600, 500))
end

function sort_spectrum(λ)
    idx_pos = imag.(λ) .> 0
    pos = λ[idx_pos]
    neg = λ[.!(idx_pos)]
    sort!(pos; lt = !isless, by = real)
    sort!(neg; lt = isless, by = real)
    return vcat(pos, neg)
end

function compute_spectrum_1d(; initial_refinement_level, nnodes,
                               accuracy_order)
    equations = LinearScalarAdvectionEquation1D(1.0)

    function initial_condition(x, t, equations::LinearScalarAdvectionEquation1D)
        return SVector(sinpi(x[1] - equations.advection_velocity[1] * t))
    end

    D_upw = upwind_operators(
        SummationByPartsOperators.Mattsson2017;
        derivative_order = 1,
        accuracy_order,
        xmin = -1.0, xmax = 1.0,
        N = nnodes)
    flux_splitting = splitting_lax_friedrichs
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralUpwind(flux_splitting),
                   volume_integral = VolumeIntegralUpwind(flux_splitting))

    coordinates_min = -1.0
    coordinates_max =  1.0
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max=10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
    J = jacobian_ad_forward(semi)
    λ = eigvals(J)
    return λ
end

function compute_spectrum_1d_dgsem(; initial_refinement_level, polydeg)
    equations = LinearScalarAdvectionEquation1D(1.0)

    function initial_condition(x, t, equations::LinearScalarAdvectionEquation1D)
        return SVector(sinpi(x[1] - equations.advection_velocity[1] * t))
    end

    solver = DGSEM(polydeg = polydeg, surface_flux = flux_godunov)

    coordinates_min = -1.0
    coordinates_max =  1.0
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max=10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
    J = jacobian_ad_forward(semi)
    λ = eigvals(J)
    return λ
end


################################################################################
# Local linear/energy experiments

function local_linear_stability(; latex = false)
    @info "Full upwind discretization (only D_−)"

    accuracy_orders = Int[]
    num_elements = Int[]
    num_nodes = Int[]
    max_real_part = Float64[]

    for accuracy_order in 2:7
        for nelements in 1:2
            for nnodes in 13:14
                λ = compute_spectrum_burgers_upwind_full(; accuracy_order,
                                                           nnodes,
                                                           nelements)
                push!(accuracy_orders, accuracy_order)
                push!(num_elements, nelements)
                push!(num_nodes, nnodes)
                push!(max_real_part, maximum(real, λ))
            end
        end
    end

    # print results
    data = hcat(accuracy_orders, num_elements, num_nodes, max_real_part)
    header = ["order", "#elements", "#nodes", "max. real part"]
    kwargs = (; header, formatters=(ft_printf("%2d", [1, 2, 3]),
                                    ft_printf("%9.2e", [4])))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end

    return nothing
end

function compute_spectrum_burgers_upwind_full(; accuracy_order,
                                                nnodes,
                                                nelements)

    D_local = derivative_operator(Mattsson2017(:minus);
                                  derivative_order = 1,
                                  xmin = -1.0, xmax = 1.0,
                                  accuracy_order, N = nnodes)
    mesh = UniformPeriodicMesh1D(xmin = -1.0, xmax = 1.0, Nx = nelements)

    D = couple_discontinuously(D_local, mesh, Val(:minus))

    u0 = rand(size(D, 2))
    J = Trixi.ForwardDiff.jacobian(u0) do u
        -D * (u.^2 ./ 2)
    end
    return eigvals(J)
end


################################################################################
# Convergence experiments

function convergence_tests_1d_advection(; latex = false)
    @info "1D linear advection"
    for accuracy_order in 2:5
        @show accuracy_order

        refinement_levels = 0:7
        num_nodes = fill(20, size(refinement_levels))
        _convergence_tests_1d_advection(; refinement_levels, num_nodes,
                                          accuracy_order, latex)

        num_nodes = 10 .* 2 .^ (0:7)
        refinement_levels = fill(2, size(num_nodes))
        _convergence_tests_1d_advection(; refinement_levels, num_nodes,
                                          accuracy_order, latex)
    end

    return nothing
end

function _convergence_tests_1d_advection(; refinement_levels, num_nodes, accuracy_order,
                                           latex = false)
    num_elements = Vector{Int}()
    errors = Vector{Float64}()

    for (initial_refinement_level, nnodes) in zip(refinement_levels, num_nodes)
        nelements = 2^initial_refinement_level
        tol = 1.0e-12
        res = compute_errors_1d_advection(; initial_refinement_level, nnodes,
                                            accuracy_order, tol)
        push!(num_elements, nelements)
        push!(errors, first(res.l2))
    end
    if length(unique(num_elements)) == 1
        eoc = compute_eoc(num_nodes, errors)
    elseif length(unique(num_nodes)) == 1
        eoc = compute_eoc(num_elements, errors)
    end

    # print results
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

function compute_errors_1d_advection(; initial_refinement_level, nnodes,
                                       accuracy_order, tol)
    equations = LinearScalarAdvectionEquation1D(1.0)

    function initial_condition(x, t, equations::LinearScalarAdvectionEquation1D)
        return SVector(sinpi(x[1] - equations.advection_velocity[1] * t))
    end

    D_upw = upwind_operators(
        SummationByPartsOperators.Mattsson2017;
        derivative_order = 1,
        accuracy_order,
        xmin = -1.0, xmax = 1.0,
        N = nnodes)
    flux_splitting = splitting_lax_friedrichs
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralUpwind(flux_splitting),
                   volume_integral = VolumeIntegralUpwind(flux_splitting))

    coordinates_min = -1.0
    coordinates_max =  1.0
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max=10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    ode = semidiscretize(semi, (0.0, 5.0))
    sol = solve(ode, RDPK3SpFSAL49(); ode_default_options()...,
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


function convergence_tests_1d_euler(; latex = false)
    for accuracy_order in 2:5
        refinement_levels = 0:7
        num_nodes = fill(20, size(refinement_levels))
        splitting = splitting_vanleer_haenel
        @info "1D compressible Euler equations" accuracy_order splitting
        _convergence_tests_1d_euler(; refinement_levels, num_nodes,
                                      accuracy_order, splitting, latex)

        refinement_levels = 0:7
        num_nodes = fill(20, size(refinement_levels))
        splitting = splitting_steger_warming
        @info "1D compressible Euler equations" accuracy_order splitting
        _convergence_tests_1d_euler(; refinement_levels, num_nodes,
                                      accuracy_order, splitting, latex)

        num_nodes = 10 .* 2 .^ (0:7)
        refinement_levels = fill(2, size(num_nodes))
        splitting = splitting_steger_warming
        @info "1D compressible Euler equations" accuracy_order splitting
        _convergence_tests_1d_euler(; refinement_levels, num_nodes,
                                      accuracy_order, splitting, latex)
    end

    return nothing
end

function _convergence_tests_1d_euler(; refinement_levels, num_nodes, accuracy_order,
                                       splitting, latex = false)
    num_elements = Vector{Int}()
    errors = Vector{Float64}()

    for (initial_refinement_level, nnodes) in zip(refinement_levels, num_nodes)
        nelements = 2^initial_refinement_level
        tol = 1.0e-13
        res = compute_errors_1d_euler(; initial_refinement_level, nnodes,
                                        accuracy_order, tol, splitting)
        push!(num_elements, nelements)
        push!(errors, first(res.l2))
    end
    if length(unique(num_elements)) == 1
        eoc = compute_eoc(num_nodes, errors)
    elseif length(unique(num_nodes)) == 1
        eoc = compute_eoc(num_elements, errors)
    end

    # print results
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

function compute_errors_1d_euler(; initial_refinement_level, nnodes,
                                   accuracy_order, tol, splitting)
    equations = CompressibleEulerEquations1D(1.4)

    initial_condition = initial_condition_convergence_test
    source_terms = source_terms_convergence_test

    D_upw = upwind_operators(
        SummationByPartsOperators.Mattsson2017;
        derivative_order = 1,
        accuracy_order,
        xmin = -1.0, xmax = 1.0,
        N = nnodes)
    solver = FDSBP(D_upw,
                   surface_integral = SurfaceIntegralUpwind(splitting),
                   volume_integral = VolumeIntegralUpwind(splitting))

    coordinates_min = 0.0
    coordinates_max = 2.0
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
                                        solver; source_terms)

    ode = semidiscretize(semi, (0.0, 2.0))
    sol = solve(ode, RDPK3SpFSAL49(); ode_default_options()...,
                abstol = tol, reltol = tol)

    analysis_callback = AnalysisCallback(semi)
    return analysis_callback(sol)
end


################################################################################
# Isentropic vortex

function experiments_isentropic_vortex()
    isentropic_vortex_generate_data()
    isentropic_vortex_plot_results(version = "blocks")
    isentropic_vortex_plot_results(version = "periodic")
end

function isentropic_vortex_generate_data()
    for accuracy_order in 2:6
        @info "Generate data" accuracy_order
        t, error_density = compute_error_isentropic_vortex(; accuracy_order,
                                                             initial_refinement_level = 2,
                                                             nnodes = 16,
                                                             source_of_coefficients = Mattsson2017
                                                             )
        open(joinpath(figdir, "isentropic_vortex_order_$(accuracy_order)_blocks.dat"), "w") do io
            println(io, "# t\tL2_error_density")
            writedlm(io, hcat(t, error_density))
        end

        t, error_density = compute_error_isentropic_vortex(; accuracy_order,
                                                             initial_refinement_level = 0,
                                                             nnodes = 64,
                                                             source_of_coefficients = periodic_derivative_operator
                                                             )
        open(joinpath(figdir, "isentropic_vortex_order_$(accuracy_order)_periodic.dat"), "w") do io
            println(io, "# t\tL2_error_density")
            writedlm(io, hcat(t, error_density))
        end
    end
end

function isentropic_vortex_plot_results(; version = "blocks")
    fig = plot(xguide = L"Time $t$", yguide = L"$L^2$ error of the density";
               xscale = :log10, yscale = :log10,
               plot_kwargs()...)

    linestyles = [:solid, :dash, :dashdot, :dot, :solid]
    for(accuracy_order, linestyle) in zip(2:6, linestyles)
        data = readdlm(joinpath(figdir,  "isentropic_vortex_order_$(accuracy_order)_$(version).dat"), comments = true)
        plot!(fig, data[:, 1], data[:, 2];
              label = "Order $(accuracy_order)", linestyle,
              plot_kwargs()...)
    end

    plot!(fig, legend = :bottomright)
    savefig(fig, joinpath(figdir, "isentropic_vortex_error_$(version).pdf"))
    @info "Error plot saved in the directory `figdir`" figdir
    return nothing
end

function compute_error_isentropic_vortex(; accuracy_order = 4, nnodes = 16,
                                           initial_refinement_level = 2,
                                           flux_splitting = splitting_steger_warming,
                                           source_of_coefficients = Mattsson2017,
                                           polydeg = nothing,
                                           volume_flux = flux_ranocha_turbo,
                                           surface_flux = flux_lax_friedrichs,
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
        ϱ0 = 1.0               # background density
        v0 = SVector(1.0, 1.0) # background velocity
        p0 = 10.0              # background pressure
        ε = 10.0               # vortex strength
        L = 10.0               # size of the domain per coordinate direction

        T0 = p0 / ϱ0           # background temperature
        γ = equations.gamma    # ideal gas constant

        vortex_center(x, L) = mod(x + L/2, L) - L/2
        x0 = v0 * t            # current center of the vortex
        dx = vortex_center.(x - x0, L)
        r2 = sum(abs2, dx)

        # perturbed primitive variables
        T = T0 - (γ - 1) * ε^2 / (8 * γ * π^2) * exp(1 - r2)
        v = v0 + ε / (2 * π) * exp(0.5 * (1 - r2)) * SVector(-dx[2], dx[1])
        ϱ = ϱ0 * (T / T0)^(1 / (γ - 1))
        p = ϱ * T

        return prim2cons(SVector(ϱ, v..., p), equations)
    end

    if polydeg === nothing
        # Use upwind SBP discretization
        D_upw = upwind_operators(source_of_coefficients;
                                 derivative_order = 1,
                                 accuracy_order,
                                 xmin = -1.0, xmax = 1.0,
                                 N = nnodes)
        solver = FDSBP(D_upw,
                       surface_integral = SurfaceIntegralUpwind(flux_splitting),
                       volume_integral = VolumeIntegralUpwind(flux_splitting))
    else
        # Use DGSEM
        volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
        solver = DGSEM(; polydeg, surface_flux, volume_integral)
    end

    coordinates_min = (-5.0, -5.0)
    coordinates_max = ( 5.0,  5.0)
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 100_000,
                    periodicity = true)

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

    sol = solve(ode, RDPK3SpFSAL49();
                abstol = tol, reltol = tol,
                ode_default_options()..., callback = saving, tstops = saveat)
    return (; t = saved_values.t, error_density = saved_values.saveval)
end


################################################################################
# Kelvin-Helmholtz instability

function experiments_kelvin_helmholtz_instability(; latex = false)
    # Upwind SBP operators
    @info "Upwind SBP operators"
    nnodes = 16
    initial_refinement_levels = 0:4
    accuracy_orders = 2:7
    for flux_splitting in (splitting_vanleer_haenel, splitting_steger_warming)
        @info flux_splitting
        final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                             length(accuracy_orders))

        for (j, accuracy_order) in enumerate(accuracy_orders)
            for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
                t = blowup_kelvin_helmholtz(; accuracy_order, nnodes,
                                            initial_refinement_level,
                                            flux_splitting)
                final_times[i, j] = t
            end
        end

        # print results
        data = hcat(4 .^ initial_refinement_levels, final_times)
        header = vcat(0, accuracy_orders)
        kwargs = (; header, title = string(flux_splitting) * ", constant #nodes",
                    formatters=(ft_printf("%3d", [1]),
                                ft_printf("%.2f", 2:size(data,2))))
        pretty_table(data; kwargs...)
        if latex
            pretty_table(data; kwargs..., backend=Val(:latex))
        end
        println()
    end

    # Upwind SBP operators with constant number of degrees of freedom
    @info "Upwind SBP operators with constant number of degrees of freedom"
    initial_refinement_levels = 0:4
    accuracy_orders = 2:7
    for flux_splitting in (splitting_vanleer_haenel, splitting_steger_warming)
        @info flux_splitting
        final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                             length(accuracy_orders))

        for (j, accuracy_order) in enumerate(accuracy_orders)
            for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
                nnodes = 256 ÷ 2^initial_refinement_level
                t = blowup_kelvin_helmholtz(; accuracy_order, nnodes,
                                            initial_refinement_level,
                                            flux_splitting)
                final_times[i, j] = t
            end
        end

        # print results
        data = hcat(4 .^ initial_refinement_levels, final_times)
        header = vcat(0, accuracy_orders)
        kwargs = (; header, title = string(flux_splitting) * ", constant #DOFs",
                    formatters=(ft_printf("%3d", [1]),
                                ft_printf("%.2f", 2:size(data,2))))
        pretty_table(data; kwargs...)
        if latex
            pretty_table(data; kwargs..., backend=Val(:latex))
        end
        println()
    end

    # DGSEM
    @info "DGSEM"
    initial_refinement_levels = 2:5
    polydegs = 2:7
    for volume_flux in (flux_ranocha, flux_shima_etal)
        @info volume_flux
        final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                             length(polydegs))

        for (j, polydeg) in enumerate(polydegs)
            for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
                t = blowup_kelvin_helmholtz(; polydeg,
                                              initial_refinement_level,
                                              volume_flux)
                final_times[i, j] = t
            end
        end

        # print results
        data = hcat(4 .^ initial_refinement_levels, final_times)
        header = vcat(0, polydegs)
        kwargs = (; header, title = string(volume_flux),
                    formatters=(ft_printf("%3d", [1]),
                                ft_printf("%.2f", 2:size(data,2))))
        pretty_table(data; kwargs...)
        if latex
            pretty_table(data; kwargs..., backend=Val(:latex))
        end
        println()
    end

    # Periodic upwind SBP operators
    @info "Periodic upwind SBP operators"
    numbers_of_nodes = 16 .* 2 .^ (0:4)
    initial_refinement_level = 0
    accuracy_orders = 2:7
    source_of_coefficients = periodic_derivative_operator
    for flux_splitting in (splitting_vanleer_haenel, splitting_steger_warming)
        @info flux_splitting source_of_coefficients
        final_times = Matrix{Float64}(undef, length(numbers_of_nodes),
                                             length(accuracy_orders))

        for (j, accuracy_order) in enumerate(accuracy_orders)
            for (i, nnodes) in enumerate(numbers_of_nodes)
                t = blowup_kelvin_helmholtz(; accuracy_order, nnodes,
                                            initial_refinement_level,
                                            source_of_coefficients,
                                            flux_splitting)
                final_times[i, j] = t
            end
        end

        # print results
        data = hcat(numbers_of_nodes .^ 2, final_times)
        header = vcat(0, accuracy_orders)
        kwargs = (; header, title = string(flux_splitting) * ", periodic upwind SBP",
                    formatters=(ft_printf("%3d", [1]),
                                ft_printf("%.2f", 2:size(data,2))))
        pretty_table(data; kwargs...)
        if latex
            pretty_table(data; kwargs..., backend=Val(:latex))
        end
        println()
    end

    return nothing
end

function blowup_kelvin_helmholtz(; accuracy_order = 4, nnodes = 16,
                                   initial_refinement_level = 2,
                                   flux_splitting = splitting_vanleer_haenel,
                                   source_of_coefficients = Mattsson2017,
                                   polydeg = nothing,
                                   volume_flux = flux_ranocha,
                                   tol = 1.0e-6)
    equations = CompressibleEulerEquations2D(1.4)

    function initial_condition(x, t, equations::CompressibleEulerEquations2D)
        # change discontinuity to tanh
        # typical resolution 128^2, 256^2
        # domain size is [-1,+1]^2
        slope = 15
        B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
        rho = 0.5 + 0.75 * B
        v1 = 0.5 * (B - 1)
        v2 = 0.1 * sin(2 * pi * x[1])
        p = 1.0
        return prim2cons(SVector(rho, v1, v2, p), equations)
    end

    if polydeg === nothing
        # Use upwind SBP discretization
        D_upw = upwind_operators(source_of_coefficients;
                                 derivative_order = 1,
                                 accuracy_order,
                                 xmin = -1.0, xmax = 1.0,
                                 N = nnodes)
        solver = FDSBP(D_upw,
                       surface_integral = SurfaceIntegralUpwind(flux_splitting),
                       volume_integral = VolumeIntegralUpwind(flux_splitting))

        @info "Kelvin-Helmholtz instability" accuracy_order nnodes initial_refinement_level flux_splitting
    else
        # Use DGSEM
        volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
        solver = DGSEM(; polydeg, surface_flux = flux_lax_friedrichs, volume_integral)

        @info "Kelvin-Helmholtz instability" polydeg initial_refinement_level volume_flux
    end

    coordinates_min = (-1.0, -1.0)
    coordinates_max = ( 1.0,  1.0)
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 100_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
    @show Trixi.ndofs(semi)

    tspan = (0.0, 15.0)
    ode = semidiscretize(semi, tspan)

    integrator = init(ode, SSPRK43(); controller = PIDController(0.55, -0.27, 0.05),
                      abstol = tol, reltol = tol,
                      ode_default_options()...)

    try
        solve!(integrator)
    catch error
        @info "Blow-up" integrator.t
        reset_threads!()
    end

    return integrator.t
end

# https://github.com/JuliaSIMD/Polyester.jl/issues/30
function reset_threads!()
    PolyesterWeave.reset_workers!()
    for i in 1:(Threads.nthreads() - 1)
        ThreadingUtilities.initialize_task(i)
    end
    return nothing
end

function run_kelvin_helmholtz(; accuracy_order = 4, nnodes = 16,
                                initial_refinement_level = 2,
                                flux_splitting = splitting_vanleer_haenel,
                                source_of_coefficients = Mattsson2017,
                                polydeg = nothing,
                                volume_flux = flux_ranocha_turbo,
                                tol = 1.0e-6)
    equations = CompressibleEulerEquations2D(1.4)

    function initial_condition(x, t, equations::CompressibleEulerEquations2D)
        # change discontinuity to tanh
        # typical resolution 128^2, 256^2
        # domain size is [-1,+1]^2
        slope = 15
        B = tanh(slope * x[2] + 7.5) - tanh(slope * x[2] - 7.5)
        rho = 0.5 + 0.75 * B
        v1 = 0.5 * (B - 1)
        v2 = 0.1 * sin(2 * pi * x[1])
        p = 1.0
        return prim2cons(SVector(rho, v1, v2, p), equations)
    end

    if polydeg === nothing
        # Use upwind SBP discretization
        D_upw = upwind_operators(source_of_coefficients;
                                 derivative_order = 1,
                                 accuracy_order,
                                 xmin = -1.0, xmax = 1.0,
                                 N = nnodes)
        solver = FDSBP(D_upw,
                       surface_integral=SurfaceIntegralUpwind(flux_splitting),
                    #    surface_integral=SurfaceIntegralStrongForm(flux_lax_friedrichs),
                       volume_integral=VolumeIntegralUpwind(flux_splitting))
    else
        # Use DGSEM
        volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
        solver = DGSEM(; polydeg, surface_flux = flux_lax_friedrichs, volume_integral)
    end

    coordinates_min = (-1.0, -1.0)
    coordinates_max = ( 1.0,  1.0)
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 100_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    tspan = (0.0, 15.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 1000
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    saving_callback = SaveSolutionCallback(; interval = 100,
                                             save_final_solution = true,
                                             output_directory = joinpath(@__DIR__, "out_dev"),
                                            #  solution_variables = cons2cons)
                                             solution_variables = cons2prim)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback,
                            alive_callback,
                            saving_callback)

    integrator = init(ode, SSPRK43(); controller = PIDController(0.55, -0.27, 0.05),
                      abstol = tol, reltol = tol,
                      ode_default_options()..., callback=callbacks)
    try
        solve!(integrator)
    catch err
        @warn "Crashed at time" integrator.t
        saving_callback.affect!(integrator)
        reset_threads!()
    end
    summary_callback() # print the timer summary

    return nothing
end


################################################################################
# Inviscid Taylor-Green vortex

function blowup_experiments_taylor_green_vortex(; latex = false)
    # Upwind SBP operators with constant number of nodes per element
    @info "Upwind SBP operators with constant number of nodes per element"
    nnodes = 16
    initial_refinement_levels = 0:2
    accuracy_orders = 2:7
    flux_splitting = splitting_steger_warming
    @info flux_splitting
    final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                         length(accuracy_orders))

    for (j, accuracy_order) in enumerate(accuracy_orders)
        for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
            t = blowup_taylor_green(; accuracy_order, nnodes,
                                      initial_refinement_level,
                                      flux_splitting)
            final_times[i, j] = t
        end
    end

    # print results
    data = hcat(8 .^ initial_refinement_levels, final_times)
    header = vcat(0, accuracy_orders)
    kwargs = (; header, title = string(flux_splitting) * ", constant #nodes",
                formatters=(ft_printf("%3d", [1]),
                            ft_printf("%.2f", 2:size(data,2))))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    println()

    # Upwind SBP operators with constant number of degrees of freedom
    @info "Upwind SBP operators with constant number of degrees of freedom"
    initial_refinement_levels = 0:2
    accuracy_orders = 2:7
    flux_splitting = splitting_steger_warming
    @info flux_splitting
    final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                        length(accuracy_orders))

    for (j, accuracy_order) in enumerate(accuracy_orders)
        for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
            nnodes = 64 ÷ 2^initial_refinement_level
            t = blowup_taylor_green(; accuracy_order, nnodes,
                                      initial_refinement_level,
                                      flux_splitting)
            final_times[i, j] = t
        end
    end

    # print results
    data = hcat(8 .^ initial_refinement_levels, final_times)
    header = vcat(0, accuracy_orders)
    kwargs = (; header, title = string(flux_splitting) * ", constant #DOFs",
                formatters=(ft_printf("%3d", [1]),
                            ft_printf("%.2f", 2:size(data,2))))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    println()

    # Periodic upwind SBP operators
    @info "Periodic upwind SBP operators"
    numbers_of_nodes = 16 .* 2 .^ (0:2)
    initial_refinement_level = 0
    accuracy_orders = 2:7
    source_of_coefficients = periodic_derivative_operator
    flux_splitting = splitting_steger_warming
    @info flux_splitting source_of_coefficients
    final_times = Matrix{Float64}(undef, length(numbers_of_nodes),
                                         length(accuracy_orders))

    for (j, accuracy_order) in enumerate(accuracy_orders)
        for (i, nnodes) in enumerate(numbers_of_nodes)
            t = blowup_taylor_green(; accuracy_order, nnodes,
                                      initial_refinement_level,
                                      source_of_coefficients,
                                      flux_splitting)
            final_times[i, j] = t
        end
    end

    # print results
    data = hcat(numbers_of_nodes .^ 3, final_times)
    header = vcat(0, accuracy_orders)
    kwargs = (; header, title = string(flux_splitting) * ", periodic upwind SBP",
                formatters=(ft_printf("%3d", [1]),
                            ft_printf("%.2f", 2:size(data,2))))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    println()

    return nothing
end

function blowup_experiments_taylor_green_vortex_mach04(; latex = false)
    # Upwind SBP operators with constant number of nodes per element
    @info "Upwind SBP operators with constant number of nodes per element"
    nnodes = 16
    initial_refinement_levels = 0:2
    accuracy_orders = 2:7
    flux_splitting = splitting_steger_warming
    @info flux_splitting
    final_times = Matrix{Float64}(undef, length(initial_refinement_levels),
                                         length(accuracy_orders))

    for (j, accuracy_order) in enumerate(accuracy_orders)
        for (i, initial_refinement_level) in enumerate(initial_refinement_levels)
            t = blowup_taylor_green(; accuracy_order, nnodes,
                                      initial_refinement_level,
                                      flux_splitting,
                                      Mach = 0.4)
            final_times[i, j] = t
        end
    end

    # print results
    data = hcat(8 .^ initial_refinement_levels, final_times)
    header = vcat(0, accuracy_orders)
    kwargs = (; header, title = string(flux_splitting) * ", constant #nodes",
                formatters=(ft_printf("%3d", [1]),
                            ft_printf("%.2f", 2:size(data,2))))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    println()

    # Periodic upwind SBP operators
    @info "Periodic upwind SBP operators"
    numbers_of_nodes = 16 .* 2 .^ (0:2)
    initial_refinement_level = 0
    accuracy_orders = 2:7
    source_of_coefficients = periodic_derivative_operator
    flux_splitting = splitting_steger_warming
    @info flux_splitting source_of_coefficients
    final_times = Matrix{Float64}(undef, length(numbers_of_nodes),
                                         length(accuracy_orders))

    for (j, accuracy_order) in enumerate(accuracy_orders)
        for (i, nnodes) in enumerate(numbers_of_nodes)
            t = blowup_taylor_green(; accuracy_order, nnodes,
                                      initial_refinement_level,
                                      source_of_coefficients,
                                      flux_splitting,
                                      Mach = 0.4)
            final_times[i, j] = t
        end
    end

    # print results
    data = hcat(numbers_of_nodes .^ 3, final_times)
    header = vcat(0, accuracy_orders)
    kwargs = (; header, title = string(flux_splitting) * ", periodic upwind SBP",
                formatters=(ft_printf("%3d", [1]),
                            ft_printf("%.2f", 2:size(data,2))))
    pretty_table(data; kwargs...)
    if latex
        pretty_table(data; kwargs..., backend=Val(:latex))
    end
    println()

    return nothing
end

function blowup_taylor_green(; accuracy_order = 4, nnodes = 16,
                               initial_refinement_level = 1,
                               flux_splitting = splitting_steger_warming,
                               source_of_coefficients = Mattsson2017,
                               polydeg = nothing,
                               volume_flux = flux_ranocha_turbo,
                               tol = 1.0e-6,
                               Mach = 0.1)
    equations = CompressibleEulerEquations3D(1.4)

    function initial_condition(x, t, equations::CompressibleEulerEquations3D)
        A  = 1.0 # magnitude of speed
        Ms = Mach # maximum Mach number

        rho = 1.0
        v1  =  A * sin(x[1]) * cos(x[2]) * cos(x[3])
        v2  = -A * cos(x[1]) * sin(x[2]) * cos(x[3])
        v3  = 0.0
        p   = (A / Ms)^2 * rho / equations.gamma # scaling to get Ms
        p   = p + 1.0/16.0 * A^2 * rho * (cos(2 * x[1]) * cos(2 * x[3]) +
                                            2 * cos(2 * x[2]) + 2 * cos(2 * x[1]) +
                                            cos(2 * x[2]) * cos(2 * x[3]))

        return prim2cons(SVector(rho, v1, v2, v3, p), equations)
    end

    if polydeg === nothing
        # Use upwind SBP discretization
        D_upw = upwind_operators(source_of_coefficients;
                                 derivative_order = 1,
                                 accuracy_order,
                                 xmin = -1.0, xmax = 1.0,
                                 N = nnodes)
        solver = FDSBP(D_upw,
                       surface_integral=SurfaceIntegralUpwind(flux_splitting),
                       volume_integral=VolumeIntegralUpwind(flux_splitting))

        @info "Taylor-Green vortex" accuracy_order nnodes initial_refinement_level flux_splitting
    else
        # Use DGSEM
        volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
        solver = DGSEM(; polydeg, surface_flux = flux_lax_friedrichs, volume_integral)

        @info "Taylor-Green vortex" polydeg initial_refinement_level volume_flux
    end

    coordinates_min = (-1.0, -1.0, -1.0) .* pi
    coordinates_max = ( 1.0,  1.0,  1.0) .* pi
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 100_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)
    @show Trixi.ndofs(semi)

    tspan = (0.0, 20.0)
    ode = semidiscretize(semi, tspan)

    integrator = init(ode, SSPRK43(); controller = PIDController(0.55, -0.27, 0.05),
                      abstol = tol, reltol = tol,
                      ode_default_options()...)

    try
        solve!(integrator)
    catch error
        @info "Blow-up" integrator.t
        reset_threads!()
    end

    return integrator.t
end

function dissipation_experiments_taylor_green_vortex()
    # Run simulations and save analysis data
    accuracy_order = 6
    initial_refinement_level = 2
    nnodes = 16
    filename = "analysis_TGV_upwind_SBP_accuracy$(accuracy_order)_level$(initial_refinement_level)_nodes$(nnodes).dat"
    run_taylor_green(; initial_refinement_level, nnodes, accuracy_order,
                       analysis_interval = 10,
                       analysis_filename = joinpath(@__DIR__, filename))

    polydeg = 3
    initial_refinement_level = 4
    filename = "analysis_TGV_DGSEM_level$(initial_refinement_level)_polydeg$(polydeg).dat"
    run_taylor_green(; initial_refinement_level, polydeg,
                       analysis_interval = 10,
                       analysis_filename = joinpath(@__DIR__, filename))

    # Plot results
    fig_kinetic_energy = plot(xguide = L"Time $t$",
                              yguide = L"Kinetic energy $E_\mathrm{kin}$")
    fig_dissipation_rate = plot(xguide = L"Time $t$",
                                yguide = L"Dissipation rate $-\Delta E_\mathrm{kin} / \Delta t$")

    kwargs = (; label = "DGSEM", plot_kwargs()...)
    data = readdlm(joinpath(@__DIR__, "analysis_TGV_DGSEM_level4_polydeg3.dat"), comments = true)
    time = data[:, 2]
    kinetic_energy = data[:, 16]
    plot!(fig_kinetic_energy, time, kinetic_energy; kwargs...)
    plot!(fig_dissipation_rate, dissipation_rate(time, kinetic_energy)...; kwargs...)

    kwargs = (; label = "Upwind SBP", linestyle = :dot, plot_kwargs()...)
    data = readdlm(joinpath(@__DIR__, "analysis_TGV_upwind_SBP_accuracy6_level2_nodes16.dat"), comments = true)
    time = data[:, 2]
    kinetic_energy = data[:, 16]
    plot!(fig_kinetic_energy, time, kinetic_energy; kwargs...)
    plot!(fig_dissipation_rate, dissipation_rate(time, kinetic_energy)...; kwargs...)

    savefig(fig_kinetic_energy, joinpath(figdir, "taylor_green_kinetic_energy.pdf"))
    savefig(fig_dissipation_rate, joinpath(figdir, "taylor_green_dissipation_rate.pdf"))

    @info "Kinetic energy plots saved in the directory `figdir`" figdir
    return nothing
end

function dissipation_rate(time, kinetic_energy)
    dissipation_rate = zeros(length(time) - 2)
    for i in eachindex(dissipation_rate)
        dissipation_rate[i] = - (kinetic_energy[i + 2] - kinetic_energy[i]) / (time[i + 2] - time[i])
    end
    return time[(begin + 1):(end - 1)], dissipation_rate
end

function run_taylor_green(; accuracy_order = 4, nnodes = 16,
                            initial_refinement_level = 1,
                            flux_splitting = splitting_steger_warming,
                            source_of_coefficients = Mattsson2017,
                            polydeg = nothing,
                            volume_flux = flux_ranocha_turbo,
                            surface_flux = flux_lax_friedrichs,
                            tol = 1.0e-6,
                            analysis_interval = 1000,
                            analysis_filename = "analysis.dat",
                            Mach = 0.1,
                            tspan = (0.0, 20.0))
    equations = CompressibleEulerEquations3D(1.4)

    function initial_condition(x, t, equations::CompressibleEulerEquations3D)
        A  = 1.0 # magnitude of speed
        Ms = Mach # maximum Mach number

        rho = 1.0
        v1  =  A * sin(x[1]) * cos(x[2]) * cos(x[3])
        v2  = -A * cos(x[1]) * sin(x[2]) * cos(x[3])
        v3  = 0.0
        p   = (A / Ms)^2 * rho / equations.gamma # scaling to get Ms
        p   = p + 1.0/16.0 * A^2 * rho * (cos(2 * x[1]) * cos(2 * x[3]) +
                                            2 * cos(2 * x[2]) + 2 * cos(2 * x[1]) +
                                            cos(2 * x[2]) * cos(2 * x[3]))

        return prim2cons(SVector(rho, v1, v2, v3, p), equations)
    end

    if polydeg === nothing
        # Use upwind SBP discretization
        D_upw = upwind_operators(source_of_coefficients;
                                 derivative_order = 1,
                                 accuracy_order,
                                 xmin = -1.0, xmax = 1.0,
                                 N = nnodes)
        solver = FDSBP(D_upw,
                       surface_integral=SurfaceIntegralUpwind(flux_splitting),
                       volume_integral=VolumeIntegralUpwind(flux_splitting))
    else
        # Use DGSEM
        volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
        solver = DGSEM(; polydeg, surface_flux, volume_integral)
    end

    coordinates_min = (-1.0, -1.0, -1.0) .* pi
    coordinates_max = ( 1.0,  1.0,  1.0) .* pi
    mesh = TreeMesh(coordinates_min, coordinates_max;
                    initial_refinement_level,
                    n_cells_max = 8 ^ (initial_refinement_level + 1))

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_callback = AnalysisCallback(semi; interval = analysis_interval,
                                         save_analysis = true,
                                         analysis_filename,
                                         extra_analysis_integrals=(energy_total,
                                                                   energy_kinetic,
                                                                   energy_internal,))

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback,
                            alive_callback)

    sol = solve(ode, SSPRK43(); controller = PIDController(0.55, -0.27, 0.05),
                abstol = tol, reltol = tol,
                ode_default_options()..., callback=callbacks)
    summary_callback() # print the timer summary
end

function performance_taylor_green_vortex()
    # number of repetitions
    n = 5

    analysis_interval = 1000
    tspan = (0.0, 1.0)

    time_rhs = Vector{Float64}(undef, n)
    for i in eachindex(time_rhs)
        run_taylor_green(; initial_refinement_level = 2, nnodes = 16,
                           accuracy_order = 6,
                           analysis_interval, tspan)
        time_rhs[i] = Trixi.TimerOutputs.time(Trixi.timer()["rhs!"]) * 1.0e-9
    end
    @info "Runtime upwind SBP rhs! in seconds" mean(time_rhs) std(time_rhs) time_rhs

    time_rhs = Vector{Float64}(undef, n)
    for i in eachindex(time_rhs)
        run_taylor_green(; initial_refinement_level = 4, polydeg = 3,
                           analysis_interval, tspan)
        time_rhs[i] = Trixi.TimerOutputs.time(Trixi.timer()["rhs!"]) * 1.0e-9
    end
    @info "Runtime DGSEM rhs! in seconds" mean(time_rhs) std(time_rhs) time_rhs

    return nothing
end
