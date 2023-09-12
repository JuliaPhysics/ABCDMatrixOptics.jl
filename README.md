# ABCDMatrixOptics

[![Build Status](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package implements the linear tracing of simple optical systems based on the [Ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis).

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
f12 = FreeSpace(200 + 300)
l2 = ThinLens(300.0)
f2 = FreeSpace(300)

# simple, single ray
beam = GeometricBeam(x=10.0, k=0.1)

# create optical system, read it from the right
M = [f2, l2, f12, l1, f1]

# apply the system to the beam
* or propagate(M, beam)
beam_p = M * beam
```





# Credits
Substantial parts of this software are forked from [ngedwin98](https://github.com/ngedwin98/) and the [ABCDBeamTrace.jl](https://github.com/ngedwin98/ABCDBeamTrace.jl/blob/master/LICENSE.md) package.
