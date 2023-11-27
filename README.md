# High-order upwind summation-by-parts methods for nonlinear conservation laws

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10200102.svg)](https://doi.org/10.5281/zenodo.10200102)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{ranocha2023high,
  title={High-order upwind summation-by-parts methods for nonlinear
         conservation laws},
  author={Ranocha, Hendrik and Winters, Andrew Ross and
          Schlottke-Lakemper, Michael and {\"O}ffner, Philipp and
          Glaubitz, Jan and Gassner, Gregor Josef},
  eprint={2311.13888},
  eprinttype={arxiv},
  eprintclass={math.NA},
  doi={10.48550/arXiv.2311.13888}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{ranocha2023highRepro,
  title={Reproducibility repository for
         "{H}igh-order upwind summation-by-parts methods for nonlinear
         conservation laws"},
  author={Ranocha, Hendrik and Winters, Andrew Ross and Schlottke-Lakemper,
          Michael and {\"O}ffner, Philipp and Glaubitz, Jan and Gassner,
          Gregor Josef},
  year={2023},
  howpublished={\url{https://github.com/trixi-framework/paper-2023-upwind}},
  doi={10.5281/zenodo.10200102}
}
```


## Abstract

High-order methods for conservation laws can be very efficient,
in particular on modern hardware. However, it can be challenging to
guarantee their stability and robustness, especially for under-resolved
flows. A typical approach is to combine a well-working baseline scheme
with additional techniques to ensure invariant domain preservation.
To obtain good results without too much dissipation, it is important to
develop suitable baseline methods.
In this article, we study upwind summation-by-parts operators, which
have been used mostly for linear problems so far. These operators come
with some built-in dissipation everywhere, not only at element interfaces
as typical in discontinuous Galerkin methods. At the same time, this
dissipation does not introduce additional parameters.
We discuss the relation of high-order upwind summation-by-parts methods
to flux vector splitting schemes and investigate their local
linear/energy stability. Finally, we present some numerical examples
for shock-free flows of the compressible Euler equations.


## Numerical experiments

The numerical experiments presented in the paper use
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl).
To reproduce the numerical experiments using Trixi.jl, you need to install
[Julia](https://julialang.org/).

The subfolder `code` of this repository contains a `README.md` file with
instructions to reproduce the numerical experiments, including postprocessing.

The numerical experiments were carried out using Julia v1.9.3.


## Authors

- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)
- [Andrew Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (RWTH Aachen University/University of Stuttgart, Germany)
- Philipp Öffner (TU Clausthal, Germany)
- Jan Glaubitz (MIT, USA)
- Gregor J. Gassner (University of Cologne, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
