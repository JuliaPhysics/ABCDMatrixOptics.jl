export AbstractBeam, GeometricBeam, GaussianBeam

"""
    abstract type AbstractBeam{T}

Serves as abstract type for `GeometricBeam` and `GaussianBeam`.
"""
abstract type AbstractBeam{T} end

"""
    GeometricBeam(w=0.0, k=0.0, z=0.0)

Defines a geometrical ray beam with keywords.
Distance to the axis is `w`, the angle with respect to the axis is `k`.
`zpos` is the initial position along the optical axis.

See also [`GaussianBeam`](@ref).

## Example
```jldoctest
julia> GeometricBeam(w=1.0)
GeometricBeam{Float64}(1.0, 0.0, 0.0)
```
"""
struct GeometricBeam{T} <: AbstractBeam{T}
    w::T
    k::T
    zpos::T
end

function GeometricBeam(; w=0.0, k=zero(w), zpos=zero(w))
    w, k, zpos = promote(w, k, zpos)
    return GeometricBeam{typeof(w)}(w, k, zpos)
end

"""
    struct GaussianBeam{T} <: AbstractBeam{T}
        q
        zpos
        n
        λ
    end
"""
struct GaussianBeam{T} <:AbstractBeam{T}
    q::Complex{T}
    zpos::T
    n::T
    λ::T
end

"""
    GaussianBeam(w0=100e-6, z=0.0, n=1.0, λ=633e-9, zpos=0.0)

Defines a physical Gaussian beam.
All parameters can be defined with keywords.

Beam initial width `w0`, propagation distance `z`, refractive index `n` 
and vacuum wavelength `λ`.
`zpos` defines a global position which is not used to calculate the q parameter but
instead is for tracking purposes

See also [`GeometricBeamBeam`](@ref).


## Example
```jldoctest
julia> GaussianBeam(w0=100e-6, z=12.0, n=1.0, λ=633e-9, zpos=0.0)
GaussianBeam{Float64}(12.0 + 0.049630215696521214im, 0.0, 1.0, 6.33e-7)
```
"""
function GaussianBeam(; w0=100e-3, z=0.0, n=1.0, λ=633e-9, zpos=0.0)
    w0, z, n, λ, zpos = promote(w0, z, n, λ, zpos)
    return GaussianBeam{typeof(w0)}(z + 1im * (π * n * w0^2) / (λ), zpos, n, λ)
end
"""
    GaussianBeam(q; zpos, n=1.0, λ=633e-9)

Returns a `GaussianBeam` defined by complex beam parameter `q`.


## Example
```jldoctest
julia> GaussianBeam(12 + 1im * 1.0 * π * 100e-6^2 / 633e-9)
GaussianBeam{Float64}(12.0 + 0.049630215696521214im, 0.0, 1.0, 6.33e-7)
```
"""
function GaussianBeam(q; zpos=0.0, n=1.0, λ=633e-9)
    zpos, n, λ = promote(zpos, n, λ)
    q = Complex{typeof(zpos)}(q)
    return GaussianBeam{real(typeof(q))}(q, zpos, n, λ)
end

"""
    z(beam::GaussianBeam)

Return the distance to the beam waist `z` of the Gaussian beam.

See also [`z`](@ref),[`zR`](@ref),[`w`](@ref),[`R`](@ref),[`w0`](@ref).
"""
z(beam::GaussianBeam) = real(beam.q)

"""
    z(beam::GaussianBeam)

Return Rayleigh length `zR` of the Gaussian beam.

See also [`z`](@ref), [`w`](@ref), [`R`](@ref), [`w0`](@ref).
"""
zR(beam::GaussianBeam) = imag(beam.q)

"""
    w0(beam::GaussianBeam)

Return `w0` of the Gaussian beam.

See also [`z`](@ref), [`zR`](@ref), [`w`](@ref), [`R`](@ref).
"""
w0(beam::GaussianBeam) = sqrt(zR(beam) * beam.λ / (π * beam.n))

"""
    w(beam::GaussianBeam)

Return current width `w` of the Gaussian beam.

See also [`z`](@ref), [`zR`](@ref), [`R`](@ref), [`w0`](@ref).
"""
w(beam::GaussianBeam) = w0(beam) * sqrt(1 + (z(beam) / zR(beam))^2) 

"""
    R(beam::GaussianBeam)

Return curvature `R` of the Gaussian beam.

See also [`z`](@ref), [`zR`](@ref), [`w`](@ref), [`w0`](@ref).
"""
R(beam::GaussianBeam) =  (z(beam) + zR(beam)^2 / z(beam)) 
