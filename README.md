# High-order upwind summation-by-parts methods for nonlinear conservation laws

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://doi.org/TODO)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{ranocha2023high,
  title={High-order upwind summation-by-parts methods for nonlinear
         conservation laws},
  author={Ranocha, Hendrik and Winters, Andrew Ross and
          Schlottke-Lakemper, Michael and {\"O}ffner, Philipp and
          Glaubitz, Jan and Gassner, Gregor Josef},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA},
  doi={TODO},
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
  howpublished={TODO},
  doi={TODO}
}
```


## Abstract

High-order methods for conservation laws can be very efficient,
in particular on modern hardware. However, it can be challenging to
guarantee their stability and robustness, especially for under-resolved
flows. A typical approach is to combine a well-working baseline scheme
with additional techniques to ensure invariant domain preservation.
To obtain good results without too much dissipation, it is important to
develop good baseline methods. Entropy-dissipative schemes using
summation-by-parts have been developed within the last decade and are
very promising --- both from a theoretical point of view and in many
applications. However, Gassner, Svärd, and Hindenlang (2022) reported
some cases where these baseline schemes fail and related this to a loss
of local linear/energy stability. In this contribution, we present
high-order upwind summation-by-parts methods as an alternative and
demonstrate their properties for under-resolved flows.


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
- Michael Schlottke-Lakemper (University of Stuttgart, Germany)
- Philipp Öffner (TU Clausthal, Germany)
- Jan Glaubitz (MIT, USA)
- Gregor J. Gassner (University of Cologne, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
