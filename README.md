<!--- Generated by README.py -->
<!--- DO NOT EDIT THIS FILE. Instead, edit and run README.py-->
# Aphrós [<img src="https://circleci.com/gh/cselab/aphros.svg?style=svg">](https://github.com/cselab/aphros/commits/master)

<img src="doc/images/foam.png" width=300 align="right">

Finite volume solver for incompressible multiphase flows with surface tension.

Key features:

- implementation in C++14
- scalability to thousands of compute nodes with the
  [Cubism](https://gitlab.ethz.ch/mavt-cse/Cubism)
  library for distributed computing on structured grids
- coroutines to enable encapsulation in the block-wise processing framework
- fluid solver based on SIMPLE or Bell-Colella-Glaz methods
- conservative split PLIC advection solver
- novel particle method for curvature estimation improving the accuracy at low resolutions
[[demo]](https://cselab.github.io/aphros/curv.html)
[[4]](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103209)

### Clone

    git clone https://github.com/cselab/aphros.git

### Documentation

<https://cselab.github.io/aphros/doc>

Generated in [doc/sphinx](doc/sphinx).

### Requirements

C++14, cmake, MPI, hdf5 (parallel), python3, python3-numpy.

Bundled dependencies:

* [hypre](https://github.com/hypre-space/hypre)
* [eigen](https://gitlab.com/libeigen/eigen) (optional)
* [overlap](https://github.com/severinstrobl/overlap) (optional)
* [vofi](https://github.com/VOFTracking/Vofi) (optional)
* [fpzip](https://github.com/LLNL/fpzip) (optional)

### Build and install

Follow instructions from [deploy/README.md](deploy/README.md) to
prepare environment and install dependencies.

Configure, build, install and run tests:

    cd src
    make -j4
    make test

### Docker

Instead of building the code in your system, you can build and run a Docker
container

    docker build github.com/cselab/aphros --tag aphros
    docker run -i aphros

### Code formatting

#### C++

Run clang-format for C/C++ files found recursively in current directory

```
ap.format
```

#### Python

Format all python files found recursively in current directory

```
yapf -ir .
```

### Add copyright notice

Add copyright notice to C/C++ source files found recursively in current
directory (if `copyright` is not found in first 10 lines of the file)

```
ap.applycopyright $(ap.findsource)
```

## Videos

Examples of simulations visualized using
[ParaView](https://www.paraview.org/) and [OSPRay](https://www.ospray.org/)
in collaboration with Jean M. Favre at [CSCS](https://www.cscs.ch).

|    |    |
:---:|:---:
[<img src="doc/images/coalescence.jpg" width=384>](https://cselab.github.io/aphros/videos/coalescence.mp4) | [<img src="doc/images/taylor_green.jpg" width=250>](https://cselab.github.io/aphros/videos/taylor_green.mp4)
Coalescence of bubbles [[4]](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103209) | Taylor-Green vortex with bubbles [[2]](https://doi.org/10.1145/3324989.3325727) [[5]](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103286)
[<img src="doc/images/vortex_bubble.jpg" width=200>](https://cselab.github.io/aphros/videos/vortex_bubble.mp4) | [<img src="doc/images/plunging_jet.jpg" width=200>](https://cselab.github.io/aphros/videos/plunging_jet.mp4)
Bubble trapped by vortex ring [[5]](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103286) | Plunging jet [[2]](https://doi.org/10.1145/3324989.3325727)
[<img src="doc/images/reactor.jpg" width=384>](https://cselab.github.io/aphros/videos/reactor.mp4) | [<img src="doc/images/mesh_bubbles.jpg" width=384>](https://cselab.github.io/aphros/videos/mesh_bubbles.mp4)
Electrochemical reactor [[1]](https://doi.org/10.1039/C9EE00219G) | Bubbles through mesh
[<img src="doc/images/rising_bubbles.jpg" width=384>](https://cselab.github.io/aphros/videos/rising_bubbles.mp4) | [<img src="doc/images/foaming_waterfall.jpg" width=384>](https://cselab.github.io/aphros/videos/foaming_waterfall.mp4)
 Rising bubbles clustering on the surface [[6]](https://doi.org/10.1103/APS.DFD.2019.GFM.V0018) [[7]](https://www.cscs.ch/publications/annual-reports/cscs-annual-report-2019) | Foaming waterfall [[8]](https://doi.org/10.1145/3394277.3401856)

|     |
|:---:|
|[<img src="doc/images/breaking_waves.jpg" width=795>](https://www.youtube.com/watch?v=iGdphpztCJQ)|
|APS Gallery of Fluid Motion 2019 award winner: Breaking waves: to foam or not to foam? [[6]](https://doi.org/10.1103/APS.DFD.2019.GFM.V0018)|

## Developers

Aphros is developed and maintained by researchers at ETH Zurich

* [Petr Karnakov](https://www.cse-lab.ethz.ch/member/petr-karnakov/)
* [Dr. Sergey Litvinov](https://www.cse-lab.ethz.ch/member/sergey-litvinov/)
* [Fabian Wermelinger](https://www.cse-lab.ethz.ch/member/fabian-wermelinger/)

under the supervision of

* [Prof. Petros Koumoutsakos](https://www.cse-lab.ethz.ch/member/petros-koumoutsakos/)

## Publications

1. S. M. H. Hashemi, P. Karnakov, P. Hadikhani, E. Chinello, S.
  Litvinov, C.  Moser, P. Koumoutsakos, and D. Psaltis, "A versatile and
  membrane-less electrochemical reactor for the electrolysis of water and
  brine", _Energy & environmental science_, 2019
  [10.1039/C9EE00219G](https://doi.org/10.1039/C9EE00219G)
2. P. Karnakov, F. Wermelinger, M. Chatzimanolakis, S. Litvinov,
  and P.  Koumoutsakos, "A high performance computing framework for multiphase,
  turbulent flows on structured grids" in _Proceedings of the platform for
  advanced scientific computing conference on – PASC ’19_, 2019
  [10.1145/3324989.3325727](https://doi.org/10.1145/3324989.3325727)
3. P. Karnakov, S. Litvinov, P. Koumoutsakos
  "Coalescence and transport of bubbles and drops"
  _10th International Conference on Multiphase Flow (ICMF)_, 2019
4. P. Karnakov, S. Litvinov, and P. Koumoutsakos, "A hybrid
  particle volume-of-fluid method for curvature estimation in multiphase
  flows”, _International journal of multiphase flow_, 2020
  [10.1016/j.ijmultiphaseflow.2020.103209](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103209)
5. Z. Wan, P. Karnakov, P. Koumoutsakos, T. Sapsis, "Bubbles in
  Turbulent Flows: Data-driven, kinematic models with history terms”,
  _International journal of multiphase flow_, 2020
  [10.1016/j.ijmultiphaseflow.2020.103286](https://doi.org/10.1016/j.ijmultiphaseflow.2020.103286)
6. P. Karnakov, S. Litvinov, J. M. Favre, P. Koumoutsakos
  "V0018: Breaking waves: to foam or not to foam?"
  _Gallery of Fluid Motion Award_
  [10.1103/APS.DFD.2019.GFM.V0018](https://doi.org/10.1103/APS.DFD.2019.GFM.V0018)
7. Annual report 2019 of the Swiss National Supercomputing Centre (cover page)
  [[link]](https://www.cscs.ch/publications/annual-reports/cscs-annual-report-2019)
8. P. Karnakov, F. Wermelinger, S. Litvinov,
  and P.  Koumoutsakos, "Aphros: High Performance Software for Multiphase Flows with Large Scale
  Bubble and Drop Clusters" in _Proceedings of the platform for
  advanced scientific computing conference on – PASC ’20_, 2020
  [10.1145/3394277.3401856](https://doi.org/10.1145/3394277.3401856)
