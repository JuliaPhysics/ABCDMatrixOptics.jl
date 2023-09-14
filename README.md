# ABCDMatrixOptics

# Early Stage, in a couple of weeks it should be stable


[![Build Status](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/JuliaPhysics/ABCDMatrixOptics.jl/graph/badge.svg?token=BHHxKcucdi)](https://codecov.io/gh/JuliaPhysics/ABCDMatrixOptics.jl)  [![Documentation for stable version](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPhysics.github.io/ABCDMatrixOptics.jl/stable) [![Documentation for development version](https://img.shields.io/badge/docs-main-blue.svg)](https://JuliaPhysics.github.io/ABCDMatrixOptics.jl/dev)

This package implements the linear tracing of simple optical systems based on the [Ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis).
The convention of the optical elements is identical to those on Wikipedia.

## Installation
Not registered yet, so
```julia
julia> ]add https://github.com/JuliaPhysics/ABCDMatrixOptics.jl
```

## Simple Example

A simple 4f re-imaging system with a magnification of 2 can be expressed as:
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
```


## To-Dos
* [x] Plotting mechanism
* [ ] Handle efractive index inside `GaussianBeam` properly (check that it works)
* [ ] Host examples as nice Pluto examples.
* [ ] Mirrors
* [ ] Plotting of Optical Elemeng such as lensesss
* [ ] Proper docstrings for all functions
* [ ] More testing
* [ ] Sagital and tangential elements
* [ ] Examples on cavities and mode matching

## Contributing
For any discussion and issues, please open an issue here. If your proposed changes are small, you can directly create a PR.
If the proposed changes are larger, it might save time on all sides if we discuss the issue beforehand!

## Credits
Substantial parts of this software are forked from [ngedwin98](https://github.com/ngedwin98/) and the [ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl) package.
