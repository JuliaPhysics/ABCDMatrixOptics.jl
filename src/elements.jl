export Elements
export Element, FreeSpace, Interface, ThinLens, ThickLens, Mirror
export Element3d, Interface3d, ThinLens3d, ThickLens3d, Mirror3d
export transfer_matrix

abstract type Elements{T} end
abstract type Element{T} <: Elements{T} end
abstract type Element3d{T} <: Elements{T} end


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
struct ThickLens{T<:Number} <: Element{T}
    R1::T
    R2::T
    t::T
    n_lens::T
    n1::T
    n2::T
end
    
function ThickLens(; R1, R2, t, n_lens=1.5, n1=1.0, n2=1.0)
    R1, R2, t, n_lens, n1, n2 = promote(R1, R2, t, n_lens, n1, n2) 
    return ThickLens{typeof(R1)}(R1, R2, t, n_lens, n1, n2)
end

"""
    ThickLens3d(θ, R1, R2, t; n_lens=1.5 n1=1.0, n2=1.0)

Construct a `θ`-tilted thick lens with the keywords:
* `R1` radius of curvature of first surface
* `R2` radius of curvature of second surface
* `t`: thickness of lens
* `n_lens=1.5` refractive index of lens
* `n1=1`: refractive index of the medium of the incoming side
* `n2=1`: refractive index of the medium of the exiting side

"""
struct ThickLens3d{T<:Number} <: Element3d{T}
	θ::T
	R1::T
    R2::T
    t::T
    n_lens::T
    n1::T
    n2::T
end
function ThickLens3d(θ; R1, R2, t, n_lens=1.5, n1=1.0, n2=1.0)
    θ, R1, R2, t, n_lens, n1, n2 = promote(θ, R1, R2, t, n_lens, n1, n2) 
    return ThickLens3d{typeof(R1)}(θ, R1, R2, t, n_lens, n1, n2)
end

struct Interface{T<:Number} <: Element{T}
    n1::T
    n2::T
    R::T
end

"""
    Interface(n1, n2)

Creates a flat interface with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(n1::Integer, n2::Integer) = Interface{Float64}(promote(n1, n2, Inf)...)
Interface(n1, n2, R=Inf) = Interface{promote_type(typeof(n1), typeof(n2), typeof(R))}(promote(n1, n2, R)...)

"""
    Interface(; n1, n2, R=Inf)

Creates a curved interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium.
"""
Interface(; n1, n2, R=Inf) = Interface{promote_type(typeof(n1), typeof(n2), typeof(R))}(promote(n1, n2, R)...)

struct Interface3d{T<:Number} <: Element3d{T}
	θ::T
    n1::T
    n2::T
    R::T
end

"""
    Interface3d(θ, n1=1, n2=1.5, R=Inf)

Creates a θ (rad.) tilted interface with radius `R` and with refractive index `n1` on the entering side
and `n2` on the new medium. `R`=Inf leads to a flat interface.
"""
Interface3d(θ ; n1=1, n2=1.5, R=Inf) = Interface3d{promote_type(typeof(θ), typeof(n1), typeof(n2), typeof(R))}(promote(θ, n1, n2, R)...)


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
	ThinLens3d(θ, f)

Creates a `θ`-tilted thin lens with focal length `f`.
"""
struct ThinLens3d{T<:Number} <: Element3d{T}
    θ::T
	f::T
end

"""
    ThinLens3d(θ; R1, R2, n_lens = 1.5, n = 1)

Creates a `θ`-tilted thin lens defined by the first radius of curvature `R1`, the second `R2`.
The lens refractive index is `n_lens` and the outer refractive index is `n`.
R>0 for concave surface
"""
function ThinLens3d(θ; R1, R2, n_lens=1.5, n=1.0) 
	f = inv((n_lens - n) / n * (-1/R1 + 1/R2))
	return ThinLens3d(θ, f)
end
ThinLens3d(θ, f) = ThinLens3d(promote(θ, f)...)

"""
    Mirror(R=Inf)

Mirror with radius of curvature `R`.
Per default `Inf`, so a flat mirror.
"""
@with_kw_noshow struct Mirror{T} <: Element{T}
    R::T=Inf
end

"""
    Mirror3d(θ, R=Inf)

θ-tilted mirror with radius of curvature `R`.
Per default `Inf`, so a flat mirror.
"""
@with_kw_noshow struct Mirror3d{T} <: Element{T}
    θ::T
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

dz(e::Interface3d{T}) where T = zero(T)
dz(e::ThinLens3d{T}) where T = zero(T)
dz(e::Mirror3d{T}) where T = zero(T)
dz(e::ThickLens3d) = e.t



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
    transfer_matrix(element::Element3d)
	
Returns both the Ray Transfer matrices (ABCDsag, ABCDtag) associated with the given, optical element
"""
# https://doi.org/10.1364/AO.26.000427 R>0 for concave surface, with e.θ2 asin(e.n1 * sin(e.θ) / e.n2) (Snell-Descartes law)
transfer_matrix(e::Interface3d) = ([1 0 ; (e.n2*cos(asin(e.n1 * sin(e.θ) / e.n2)) - e.n1*cos(e.θ) )/(e.n2*e.R) e.n1/e.n2 ],
				[cos(asin(e.n1 * sin(e.θ) / e.n2))/cos(e.θ) 0 ; ((e.n2*cos(asin(e.n1 * sin(e.θ) / e.n2)) - e.n1*cos(e.θ)) / (e.R * e.n2))  (e.n1*cos(e.θ))/(e.n2*cos(asin(e.n1 * sin(e.θ) / e.n2)))] )
transfer_matrix(e::ThinLens3d) = ([1 0 ; -1/(e.f * cos(e.θ)) 1],
				[1 0 ; -cos(e.θ)/e.f 1])
transfer_matrix(e::Mirror3d) = ([1 0 ; -2*cos(e.θ)/e.R 1] ,
				[1 0 ; -2/(e.R * cos(e.θ)) 1] )
transfer_matrix(e::ThickLens3d) = transfer_matrix([Interface3d(e.θ, n1=e.n1, n2 = e.n_lens, R=e.R1),
							FreeSpace(e.t),
							Interface3d(e.θ, n1=e.n_lens, n2 = e.n2, R=e.R2)])

"""
    transfer_matrix(elements)

Returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.
"""
function transfer_matrix(elements::Vector{<:Element})
    return mapfoldr(transfer_matrix, (a,b) -> b * a, elements)
end
function transfer_matrix(elements::Vector{<:Elements})
    return mapfoldr(transfer_matrix, (a,b) -> transfer_2matrices_3d(a,b), elements)
	return
end

function transfer_2matrices_3d(A::Tuple, B::Tuple)
	# A : (ABCDsag, ABCDtan) , B : (ABCDsag, ABCDtan)
	return (B[1] * A[1] , B[2]* A[2] )
end
function transfer_2matrices_3d(A::Matrix, B::Tuple)
	# A : ABCD , B : (ABCDsag, ABCDtan)
	return (B[1] * A , B[2] * A )
end
function transfer_2matrices_3d(A::Tuple, B::Matrix)
	# A : (ABCDsag, ABCDtan) , B : ABCD 
	return (B * A[1] , B * A[2])
end
function transfer_2matrices_3d(A::Matrix, B::Matrix)
	#Both are single ABCD matrix
	return (B * A, B * A)
end

"""
    discretize(e)

Discretizes the elements for plots. Nothing is done expect for FreeSpace, which is split up
"""

discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.dz/N), N)
discretize(e::Elements, N::Int) = e
discretize(els::Vector{<:Elements}, N::Int) = vcat(discretize.(els, Ref(N))...)


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
