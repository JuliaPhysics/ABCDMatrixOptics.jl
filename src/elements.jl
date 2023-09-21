export Element, FreeSpace, Interface, ThinLens, ThickLens, Mirror
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
    ThickLens(;  R1, R2, t, n_lens=1.5 n1=1.0, n2=1.0)

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
    
function ThickLens(; R1, R2, t, n_lens=1.5, n1=1.0, n2=1.0)
    R1, R2, t, n_lens, n1, n2 = promote(R1, R2, t, n_lens, n1, n2) 
    return ThickLens{typeof(R1)}(R1, R2, t, n_lens, n1, n2)
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
Interface(n1::Integer, n2::Integer) = Interface{Float64}(promote(n1, n2, 0.0, Inf)...)
Interface(n1::Integer, n2::Integer, R) = Interface{Float64}(promote(n1, n2, 0.0, R)...)
Interface(n1, n2) = Interface{typeof(n1)}(n1, n2, 0.0, Inf)

"""
    Interface(n1, n2, R)

Creates a curved interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1, n2, R) = Interface{typeof(n1)}(n1, n2, 0, R)


@with_kw_noshow struct ThinLens{T<:Number} <: Element{T}
    f::T
    θ::T=typeof(f)(0)
end

"""
    ThinLens(f)

Creates a thin lens with focal length `f`.
"""
ThinLens(f::T) where T = ThinLens{T}(f, 0)

"""
    ThinLens(R1, R2, n_lens, n)

Creates a thin lens defined by the first radius of curvature `R1`, the second `R2`.
The lens refractive index is `n_lens` and the outer refractive index is `n`.
"""
ThinLens(R1, R2, n_lens=1.5, n=1.0) = ThinLens(inv((n_lens - n) / n * (1/R1 - 1/R2)))

"""
    Mirror(R=Inf)

Mirror with radius of curvature `R`.
Per default `Inf`, so a flat mirror.
"""
@with_kw_noshow struct Mirror{T} <: Element{T}
    R::T=Inf
end


# definitions of dz
"""
    dz(element::Element)

Returns how much an element changes the optical distance `z`.
"""
dz(e::FreeSpace) = e.dz
dz(e::Interface{T}) where T = zero(T)
dz(e::ThinLens{T}) where T = zero(T)
dz(e::Mirror{T}) where T = zero(T)
dz(e::ThickLens) = e.t
dz(e::Matrix) = Inf




"""
    transfer_matrix(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
transfer_matrix(e::Matrix) = e
transfer_matrix(e::Interface) = [1 0 ; ((e.n1 - e.n2) / (e.R * e.n2))  (e.n1 / e.n2)]
transfer_matrix(e::ThinLens) = [1 0 ; -1/e.f 1]
transfer_matrix(e::Mirror) = [1 0 ; -2/e.R 1]
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


"""
    discretize(e)

Discretizes the elements for plots. Nothing is done expect for FreeSpace, which is split up
"""

discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.dz/N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(els, Ref(N))...)


"""
    Base.isapprox(a::Vector{<:Element}, b::Vector{<:Element})

Compare two vectors of elements using Base.isapprox for each element's
ray matrix (ABCD entries). Does consequently not consider one
discretization of element FreeSpace different from another, or one
realization of an imaging system from another as long as both achieve
(within tolerances) the same imaging.

!!! note
    The `atol` (absolute tolerance) parameter can be used but is
    typically nonsensical as it will be used for each of the
    ray matrix entries ABCD which usually differ vastly in magnitude.

"""
Base.isapprox(
              a::Union{Element ,Vector{<:Element}}, b::Union{Element, Vector{<:Element}}; kwargs...
) = isapprox(transfer_matrix(a), transfer_matrix(b); kwargs...)

Base.isapprox(a::Matrix, b::Union{Element, Vector{<:Element}}; kwargs...) = isapprox(a, transfer_matrix(b); kwargs...)
Base.isapprox(a::Union{Element, Vector{<:Element}}, b::Matrix; kwargs...) = Base.isapprox(b, a) 
