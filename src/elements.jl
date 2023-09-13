export Element, FreeSpace, Interface, ThinLens
export transfer_matrix

abstract type Element{T} end



@with_kw_noshow struct FreeSpace{T<:Number} <: Element{T}
    dz::T
end


@with_kw_noshow struct Interface{T<:Number} <: Element{T}
    n1::T=1.0
    n2::T
    θ::T=zero(n1)
    R::T=typeof(n1)(Inf) # R > 0 when light hits concave side
end

"""
    Interface(n1, n2)

Creates a flat interface with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1::T, n2::T) where T = Interface{T}(n1, n2, zero(T), T(Inf))
Interface(n1::Integer, n2::Integer) where T = Interface{Float64}(n1, n2, zero(T), T(Inf))
Interface(n1, n2) = Interface{Float64}(n1, n2, 0.0, Inf)

"""
    Interface(n1, n2, R)

Creates a curved interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1, n2, R) = Interface{Float64}(n1, n2, θ, R)


@with_kw_noshow struct ThinLens{T<:Number} <: Element{T}
    f::T
    θ::T
end

"""
    ThinLens(f)

Creates a thin lens with focal length `f`.
"""
ThinLens(f::T) where T = ThinLens{T}(f,0)


# definitions of dz
"""
    dz(element::Element)

Returns how much an element changes the optical distance `z`.
"""
dz(e::FreeSpace) = e.dz
dz(e::Interface{T}) where T = zero(T)
dz(e::ThinLens{T}) where T = zero(T)




"""
    transfer_matrix(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
transfer_matrix(e::Interface) = [1 0 ; ((e.n1 - e.n2) / (e.R * e.n2))  (e.n1 / e.n2)]
transfer_matrix(e::ThinLens, nprev=1) = [1 0 ; -1/e.f 1]
transfer_matrix(e::FreeSpace, nprev=1) = [1 e.dz ; 0 1]



"""
    transfer_matrix(elements)

Returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.
"""
function transfer_matrix(elements::Vector{<:Element})
    return mapfoldr(transfer_matrix, (a,b) -> b * a, elements)
end
