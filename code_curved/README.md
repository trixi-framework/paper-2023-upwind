# Curvilinear numerical experiments

This directory contains all source code and mesh files required to reproduce
the numerical experiments for the curvilinear upwind SBP solver presented
in the paper. It is developed for Julia v1.10.0.

To reproduce the numerical experiments, start Julia in this directory and
execute the following commands in the Julia REPL. We recommend starting
Julia with several threads to improve the runtime performance.

```julia
julia> include("code.jl")

julia> convergence_tests_curvilinear_2d_euler()

julia> free_stream_preservation_generate_data()

julia> isentropic_vortex_generate_data()
```

The solution plots of the isentropic Euler vortex are obtained by
running the following commands and visualizing the resulting
`*.vtu` files using Paraview. To reproduce the figures presented
one should use the ParaView state file `warped_vortex.pvsm`.

```julia
julia> include("code.jl")

julia> # one run - create figures after this block
       begin
       rm("out_dev", recursive = true, force = true)
       run_isentropic_vortex()
       trixi2vtk("out_dev/solution_*.h5"; output_directory = "out_dev",
         reinterpolate = false, data_is_uniform = true)
       end
```
