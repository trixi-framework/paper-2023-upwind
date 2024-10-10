# On the robustness of high-order upwind summation-by-parts methods for nonlinear conservation laws

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10200102.svg)](https://doi.org/10.5281/zenodo.10200102)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@article{ranocha2025robustness,
  title={On the robustness of high-order upwind summation-by-parts
         methods for nonlinear conservation laws},
  author={Ranocha, Hendrik and Winters, Andrew Ross and
          Schlottke-Lakemper, Michael and {\"O}ffner, Philipp and
          Glaubitz, Jan and Gassner, Gregor Josef},
  journal={Journal of Computational Physics},
  volume={520},
  pages={113471},
  year={2025},
  month={01},
  doi={10.1016/j.jcp.2024.113471},
  eprint={2311.13888},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{ranocha2025robustnessRepro,
  title={Reproducibility repository for
         "{O}n the robustness of high-order upwind summation-by-parts methods
         for nonlinear conservation laws"},
  author={Ranocha, Hendrik and Winters, Andrew Ross and Schlottke-Lakemper,
          Michael and {\"O}ffner, Philipp and Glaubitz, Jan and Gassner,
          Gregor Josef},
  year={2023},
  howpublished={\url{https://github.com/trixi-framework/paper-2023-upwind}},
  doi={10.5281/zenodo.10200102}
}
```


## Abstract

We use the framework of upwind summation-by-parts (SBP) operators developed
by Mattsson (2017, [DOI: 10.1016/j.jcp.2017.01.042](https://doi.org/10.1016/j.jcp.2017.01.042))
study different flux vector splittings in this context. To do so, we
introduce discontinuous-Galerkin-like interface terms for multi-block upwind
SBP methods applied to nonlinear conservation laws.
We investigate the behavior of the upwind SBP methods for flux vectors splittings
of varying complexity on Cartesian as well as unstructured curvilinear multi-block meshes.
Moreover, we analyze the local linear/energy stability of these methods following
Gassner, Svärd, and Hindenlang (2022, [DOI: 10.1007/s10915-021-01720-8](https://doi.org/10.1007/s10915-021-01720-8)).
Finally, we investigate the robustness of upwind SBP methods for challenging
examples of shock-free flows of the compressible Euler equations such as
a Kelvin-Helmholtz instability and the inviscid Taylor-Green vortex.


## Numerical experiments

The numerical experiments presented in the paper use
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
To reproduce the numerical experiments using Trixi.jl, you need to install
[Julia](https://julialang.org/).

The subfolder `code` of this repository contains a `README.md` file with
instructions to reproduce the Cartesian mesh numerical experiments and
the subfolder `code_curved` contains a `README.md` file with instructions
to reproduce the curvilinear mesh numerical experiments.
Both subfolders also include information about postprocessing.

The Cartesian mesh numerical experiments were carried out using Julia v1.9.3
and the curvilinear mesh results were carried out using Julia 1.10.0.


## Authors

- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)
- [Andrew Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (RWTH Aachen University/University of Stuttgart, Germany)
- Philipp Öffner (TU Clausthal, Germany)
- Jan Glaubitz (MIT, USA)
- Gregor J. Gassner (University of Cologne, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
