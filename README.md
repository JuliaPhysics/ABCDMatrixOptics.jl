# ABCDMatrixOptics.jl

[![Build Status](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/JuliaPhysics/ABCDMatrixOptics.jl/graph/badge.svg?token=BHHxKcucdi)](https://codecov.io/gh/JuliaPhysics/ABCDMatrixOptics.jl)  [![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPhysics.github.io/ABCDMatrixOptics.jl/stable) [![Documentation for development version](https://img.shields.io/badge/docs-main-blue.svg)](https://JuliaPhysics.github.io/ABCDMatrixOptics.jl/dev)

This package implements the linear tracing of simple optical systems based on the [Ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis).
The convention of the optical elements is identical to those on Wikipedia. So far we allow most standard optical elements and geometric and Gaussian beams.
If you need another element, feel free to open an issue!

## Installation
This package is registered, so type (literally type a `]` in the terminal to enter the package manager and leave the package manager with pressing backspace):
```julia
julia> ]add ABCDMatrixOptics.jl
```

## Simple Example
See the [Pluto.jl](https://github.com/fonsp/Pluto.jl) examples [here](examples/).
A simple 4f re-imaging system with a magnification of 2 can be expressed as.
```julia
using ABCDMatrixOptics

# create needed optical elements
f1 = FreeSpace(200)
l1 = ThinLens(200.0)
f12 = FreeSpace(200 + 400)
l2 = ThinLens(400.0)
f2 = FreeSpace(400)

# simple, single ray
beam = GeometricBeam(w=10.0, k=0.1)

# create optical system
# it's built from left to right.
M = [f1, l1, f12, l2, f2]

# apply the system to the beam
# The matrices are evaluates as f2 * l2 * f12 * l1 * f1 * beam
# * is syntactic sugar for propagate
beam_p = M * beam
GeometricBeam{Float64}(-19.999999999999996, -0.05000000000000001, 1200.0)
# beam is magnified by 2 in size


# GaussianBeam
red_beam = GaussianBeam(w0=5e-3)
blue_beam = GaussianBeam(w0=5e-3, λ=405e-9)

using Plots

plot(M, red_beam)
plot!(M, blue_beam)
```
![Simple plot](docs/assets/plot.png)

## To-Dos
* [x] Plotting mechanism
* [x] Mirrors
* [x] Proper docstrings for all functions
* [x] Handle refractive index inside `GaussianBeam` properly (check that it works)
* [x] More testing
* [ ] Host examples as nice Pluto examples.
* [ ] Plotting of Optical elements such as lenses
* [ ] Sagital and tangential elements
* [ ] Examples on cavities and mode matching

## Contributing
For any discussion and issues, please open an issue here. If your proposed changes are small, you can directly create a PR.
If the proposed changes are larger, it might save time on all sides if we discuss the issue beforehand!

## Credits
Substantial parts of this software are forked from [ngedwin98](https://github.com/ngedwin98/) and the [ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl) package.
