export AbstractBeam, GeometricBeam

abstract type AbstractBeam end

@with_kw struct GeometricBeam{T} <: AbstractBeam
    x::T # radial extent
    k = zero(x) # (angle/sin/tan of) slope of beam
    z = zero(x) # position along beam axis
    n = one(x) # index of refraction aka optical density
end
