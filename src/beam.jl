export AbstractBeam, GeometricBeam, GaussianBeam

"""
    abstract type AbstractBeam{T}

Serves as abstract type for `GeometricBeam` and `GaussianBeam`.
"""
abstract type AbstractBeam{T} end

"""
    GeometricBeam(w=0.0, k=0.0, z=0.0)

Define a geometrical ray beam.
Distance to the axis is `w`, the angle with respect to the axis is `k`.
`z` is the position along the optical axis.
"""
@with_kw_noshow struct GeometricBeam{T} <: AbstractBeam{T}
    w::T = zero(0.0) # radial extent
    k::T = zero(w) # (angle/sin/tan of) slope of beam
    z::T = zero(w) # position along beam axis
end

"""
    GaussianBeam(w=100e-6, z=0.0, n=1.0, λ=633e-9)

Defines a physical Gaussian beam.
All parameters can be defined with keywords.

Beam initial width `w`, initial propagation distance `z`, refractive index `n` 
and vacuum wavelength `λ`
"""
@with_kw_noshow struct GaussianBeam{T} <:AbstractBeam{T}
    w::T=100e-6
    z::T=typeof(w)(0.0)
    n::T=typeof(w)(1.0)
    λ::T=typeof(w)(633e-9)
end


"""
    q(beam::GaussianBeam{T}) where T

Returns the complex beam parameter `q` calculated from `beam`.
"""
q(beam::GaussianBeam{T}) where T = beam.z + 1im * zR(beam)
zR(beam::GaussianBeam{T}) where T = (π * beam.n) * beam.w^2 / beam.λ 
