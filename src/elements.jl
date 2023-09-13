export Element, FreeSpace, Interface, ThinLens, ThickLens
export transfer_matrix

abstract type Element{T} end


"""
    FreeSpace(dz)

Construct a free space propagation over distance `dz`.
"""
@with_kw_noshow struct FreeSpace{T<:Number} <: Element{T}
    dz::T
end

"""
    ThickLens(; n_lens, R1, R2, t)

Construct a thick lens with the keywords:

* `R1` radius of curvature of first surface
* `R2` radius of curvature of second surface
* `t`: thickness of lens
* `n_lens=1.5` refractive index of lens
* `n1=1`: refractive index of the medium of the incoming side
* `n2=1`: refractive index of the medium of the exiting side

"""
@with_kw_noshow struct ThickLens{T<:Number} <: Element{T}
    R1::T
    R2::T
    t::T
    n_lens::T=1.5
    n1::T=1.0
    n2::T=1.0
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
Interface(n1::Integer, n2::Integer) = Interface{Float64}(n1, n2, zero(T), T(Inf))
Interface(n1, n2) = Interface{Float64}(n1, n2, 0.0, Inf)

"""
    Interface(n1, n2, R)

Creates a curved interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1, n2, R) = Interface{Float64}(n1, n2, θ, R)


@with_kw_noshow struct ThinLens{T<:Number} <: Element{T}
    f::T
    θ::T=typeof(f)(0)
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
dz(e::ThickLens) = e.t
dz(e::Matrix) = Inf




"""
    transfer_matrix(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
transfer_matrix(e::Interface) = [1 0 ; ((e.n1 - e.n2) / (e.R * e.n2))  (e.n1 / e.n2)]
transfer_matrix(e::ThinLens) = [1 0 ; -1/e.f 1]
transfer_matrix(e::FreeSpace) = [1 e.dz ; 0 1]
transfer_matrix(e::ThickLens) = transfer_matrix([Interface(n1=e.n1, n2=e.n_lens, R=e.R1), 
                                                 FreeSpace(e.t), 
                                                 Interface(n1=e.n_lens, n2=e.n2, R=e.R2)])



"""
    transfer_matrix(elements)

Returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.
"""
function transfer_matrix(elements::Vector{<:Element})
    return mapfoldr(transfer_matrix, (a,b) -> b * a, elements)
end
