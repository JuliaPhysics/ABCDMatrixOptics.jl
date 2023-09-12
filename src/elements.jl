export Element, FreeSpace, Interface, ThinLens


abstract type Element{T} end

"""
    RTM(element::Element) 

Returns the Ray Transfer (ABCD) matrix associated with the given, optical element.
"""
RTM


@with_kw_noshow struct FreeSpace{T<:Number} <: Element{T}
    dz::T
end


RTM(e::FreeSpace, nprev=1) = [1 e.dz ; 0 1]
dz(e::FreeSpace) = e.dz


@with_kw_noshow struct Interface{T<:Number} <: Element{T}
    n::T
    θ::T=zero(n)
    R::T=typeof(x)(Inf) # R > 0 when light hits concave side
end

"""
    Interface(n)

Creates a flat interface with refractive index `n`.
"""
Interface(n::T) where T = Interface{typeof(n)}(n, zero(T), T(Inf))
Interface(n::Int) = Interface{Float64}(n, 0.0, Inf)
RTM(e::Interface, n1=1) = [1 0 ; (n1 - e.n) / e.R n1 / e.n]
dz(e::Interface{T}) where T = zero(T)


@with_kw_noshow struct ThinLens{T<:Number} <: Element{T}
    f::T
    θ::T
end

"""
    ThinLens(f)

Creates a thin lens with focal length `f`.
"""
ThinLens(f::T) where T = ThinLens{T}(f,0)
RTM(e::ThinLens, nprev=1) = [1 0 ; -1/e.f 1]
dz(e::ThinLens{T}) where T = zero(T)



"""
    RTM(elements)

Returns the Ray Transfer (ABCD) matrix associated with
an optical system described by a collection (e.g. a vector or
iteration) of optical elements.

Note, this should be used with caution. Since the output is a matrix,
it will loose information about the latest medium (refractive index) present.

Better, use a `Vector{<:Elements}` and apply it as a whole to a `b::AbstractBeam`.
"""
function RTM(elements::Vector{<:Element}, nprev = 1)
    # identity matrix
    M = [1 0; 0 1]

    # go through the elements
    # only Interface changes the refractive index permanently
    for e in reverse(elements)
        M = RTM(e, nprev) * M
        if e isa Interface 
            nprev = e.n
        end
    end
    return M
end
