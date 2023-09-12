module ABCDMatrixOptics

using Parameters
using RecipesBase

include("beam.jl")
include("elements.jl")
 #include("plots-recipes.jl")

export propagate, beamtrace




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
    a::Vector{<:Element}, b::Vector{<:Element}; kwargs...
) = isapprox(RTM(a), RTM(b); kwargs...)



"""
    propagate(e::Union{Element, Vector{<:Element}}, b)

Propagate a beam `b` either by a single element `e` or an vector
of elements.

Returned will be the final beam.

"""
function propagate(e::Element, b::GeometricBeam{T}) where T
    x, k = RTM(e) * [b.x, b.k]
    return GeometricBeam{T}(
                x=x, k=k,
                z = b.z + dz(e)
            )
end

function propagate(es::Vector{<:Element}, b::AbstractBeam)
    for e in reverse(es)
        b = propagate(e, b)
    end
    return b
end

"""
    beamtrace(elems::Vector{<:Element}, b0::AbstractBeam)

Trace a beam `b0` through a vector of elements `elems`.
All intermediate states of the beam will be recorded.

Returned will be a list of states where the last entry is the final beam.
"""
function beamtrace(elems::Vector{<:Element}, b0::B) where B
    bs = Vector{B}(undef, length(elems)+1)
    bs[1] = b0
    for (idx, elem) in enumerate(reverse(elems))
        bs[idx + 1] = propagate(elem, bs[idx])
    end
    return bs
end


end
