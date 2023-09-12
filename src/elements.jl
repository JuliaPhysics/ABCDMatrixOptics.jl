export Element, FreeSpace, Interface, ThinLens
export RTM


abstract type Element end

"""
    RTM(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
RTM


@kwdef struct FreeSpace{T<:Number} <: Element
    dz::T
end
RTM(e::FreeSpace) = [1 e.dz ; 0 1]
dz(e::FreeSpace) = e.dz


@with_kw struct Interface{T<:Number} <: Element
    n::T
    θ::T=zero(n)
    R::T=zero(n) # R > 0 when light hits concave side
end

"""
    Interface(n)

Creates a flat interface with refractive index `n`.
"""
Interface(n::T) where T = Interface{typeof(n)}(n, 0.0, Inf)
RTM(e::Interface) = [1 0 ; (e.η-1)/e.R e.η]
dz(e::Interface{T}) where T = zero(T)


struct ThinLens{T<:Number} <: Element
    f::T
    θ::T
end

"""
    ThinLens(f)

Creates a thin lens with focal length `f`.
"""
ThinLens(f::T) where T = ThinLens{T}(f,0)
RTM(e::ThinLens) = [1 0 ; -1/e.f 1]
dz(e::ThinLens{T}) where T = zero(T)



"""
    RTM(elements)

returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.
"""
RTM(elements) = mapreduce(RTM, *, elements)
