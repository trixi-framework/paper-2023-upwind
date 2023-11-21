# Numerical experiments

This directory contains all source code required to reproduce the numerical
experiments presented in the paper. It is developed for Julia v1.9.3.

To reproduce the numerical experiments, start Julia in this directory and
execute the following commands in the Julia REPL.

```julia
julia> include("code.jl")

julia> convergence_tests_1d_advection()

julia> convergence_tests_1d_euler()

julia> plot_spectra()

julia> local_linear_stability()

julia> experiments_isentropic_vortex()

julia> experiments_kelvin_helmholtz_instability()

julia> blowup_experiments_taylor_green_vortex()

julia> blowup_experiments_taylor_green_vortex_mach04()

julia> dissipation_experiments_taylor_green_vortex()

julia> performance_taylor_green_vortex()
```

The solution plots of the Kelvin-Helmholtz instability are obtained by
running one of the following commands and visualizing the resulting
`*.vtu` files using Paraview. We recommend starting Julia with several
threads to improve the runtime performance.

```julia
julia> include("code.jl")

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         accuracy_order = 6, nnodes = 16, initial_refinement_level = 4,
         flux_splitting = splitting_vanleer_haenel)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         accuracy_order = 6, nnodes = 32, initial_refinement_level = 3,
         flux_splitting = splitting_vanleer_haenel)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         accuracy_order = 6, nnodes = 64, initial_refinement_level = 2,
         flux_splitting = splitting_vanleer_haenel)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         accuracy_order = 6, nnodes = 128, initial_refinement_level = 1,
         flux_splitting = splitting_vanleer_haenel)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         accuracy_order = 6, nnodes = 256, initial_refinement_level = 0,
         flux_splitting = splitting_vanleer_haenel)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_kelvin_helmholtz(;
         polydeg = 3, initial_refinement_level = 6,
         volume_flux = flux_ranocha_turbo)
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false)
       end
```
