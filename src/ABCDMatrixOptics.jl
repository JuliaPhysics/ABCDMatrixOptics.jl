module ABCDMatrixOptics

using Parameters
using RecipesBase

include("beam.jl")
include("elements.jl")
 #include("plots-recipes.jl")

export propagate, trace




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
) = isapprox(transfer_matrix(a), transfer_matrix(b); kwargs...)



"""
    propagate(e::Union{Element, Vector{<:Element}}, b)

Propagate a beam `b` either by a single element `e` or an vector
of elements.

Returned is the final beam.

Also available as `e * b`.
"""
function propagate(e::Element, b::GeometricBeam{T}) where T
    w, k = transfer_matrix(e) * [b.w, b.k]
    return GeometricBeam{T}(w=w, k=k, z=b.z + dz(e))
end

function propagate(es::Vector{<:Element}, b::AbstractBeam)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

Base.:*(e::Union{Element, Vector{<:Element}}, b::AbstractBeam) = propagate(e, b)
Base.:*(a::Element, b::Element) = transfer_matrix(a) * transfer_matrix(b)


"""
    trace(elems::Vector{<:Element}, b0::AbstractBeam)

Trace a beam `b0` through a vector of elements `elems`.
All intermediate states of the beam will be recorded.

Return is a `Vector` of states where the last entry is the final beam.
Final beam is equivalent to `propagate(elems, b0)`.
"""
function trace(elems::Vector{<:Element}, b0::B) where B
    bs = Vector{B}(undef, length(elems)+1)
    bs[1] = b0
    for (idx, e) in enumerate(elems)
        bs[idx + 1] = propagate(e, bs[idx])
    end
    return bs
end


end
