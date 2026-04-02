# A Positivity-Preserving Relaxation Algorithm

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO.svg)](https://zenodo.org/records/10.5281/zenodo.19386973)

This repository contains information and code to reproduce the results presented
in the article
```bibtex
@online{IRS2026,
      title={A Positivity-Preserving Relaxation Algorithm}, 
      author={Thomas Izgin and Hendrik Ranocha and Chi-Wang Shu},
      year={2026},
      eprint={TODO},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      url={https://arxiv.org/abs/TODO}, 
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{IRS2026repository,
  title={Reproducibility repository for
         "A Positivity-Preserving Relaxation Algorithm"},
  author={Izgin, Thomas and Ranocha, Hendrik and Shu, Chi-Wang},
  year={2026},
  howpublished={\url{https://github.com/IzginThomas/PositiveRelaxation}},
  doi={10.5281/zenodo.19386973}
}
```

## Abstract

We combine Patankar-type methods with suitable relaxation procedures that are capable of ensuring correct dissipation or conservation of functionals such as entropy or energy while producing unconditionally positive and conservative approximations. To that end, we adapt the relaxation algorithm to enforce positivity by using either ideas from the dense output framework when a linear invariant must be preserved, or simply a geometric mean if the only constraint is positivity preservation. The latter merely requires the solution of a scalar nonlinear equation while former results in a coupled linear-nonlinear system of equations. We present sufficient conditions for the solvability of the respective equations. Several applications in the context of ordinary and partial differential equations are presented, and the theoretical findings are validated numerically.


## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Matlab](https://de.mathworks.com/products/matlab.html).
The numerical experiments presented in this article were performed using
Matlab R2025b.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Matlab in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors
- [Thomas Izgin](https://uni-kassel.de/go/izgin) (University of Kassel, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)
- [Chi-Wang Shu](https://www.dam.brown.edu/people/shu/) (Brown University, Rhode Island, USA)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
