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
    zpos::T = zero(w) # position along beam axis
end

"""
    GaussianBeam(w=100e-6, z=0.0, n=1.0, λ=633e-9)

Defines a physical Gaussian beam.
All parameters can be defined with keywords.

Beam initial width `w`, initial propagation distance `z`, refractive index `n` 
and vacuum wavelength `λ`
"""
@kwdef struct GaussianBeam{T} <:AbstractBeam{T}
    q::Complex{T}
    zpos::T=typeof(w)(0.0)
    n::T=typeof(w)(1.0)
    λ::T=typeof(w)(633e-9)
end

GaussianBeam(; w0=100e-3, zpos=0.0, n=1.0, λ=633e-9) = GaussianBeam{typeof(w0)}(1im * (π * n * w0^2) / (λ), zpos, n, λ)

"""
    q(beam::GaussianBeam{T}) where T

Returns the complex beam parameter `q` calculated from `beam`.
"""
z(beam::GaussianBeam) = real(beam.q)
zR(beam::GaussianBeam) = imag(beam.q)
w0(beam::GaussianBeam) = sqrt(zR(beam) * beam.λ / (π * beam.n))
w(beam::GaussianBeam) = w0(beam) * sqrt(1 + (z(beam) / zR(beam))^2) 
R(beam::GaussianBeam) = inv(real(inv(beam.q)))
