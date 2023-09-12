export AbstractBeam, GeometricBeam

abstract type AbstractBeam{T} end

@with_kw_noshow struct GeometricBeam{T} <: AbstractBeam{T}
    x::T = zero(0.0) # radial extent
    k::T = zero(x) # (angle/sin/tan of) slope of beam
    z::T = zero(x) # position along beam axis
    n::T = one(x) # index of refraction aka optical density
end
