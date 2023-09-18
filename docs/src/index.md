# ABCDMatrixOptics.jl

This package implements the linear tracing of simple optical systems based on the [Ray transfer matrix analysis](https://en.wikipedia.org/wiki/Ray_transfer_matrix_analysis).
The convention of the optical elements is identical to those on Wikipedia.

See on [GitHub](https://github.com/JuliaPhysics/ABCDMatrixOptics.jl) for current issues, news...

## Installation
This package is registered, so do
```julia
julia> ]add ABCDMatrixOptics.jl
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

