module ABCDMatrixOptics

using Parameters

include("beam.jl")
include("elements.jl")
include("plots-recipes.jl")

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
              a::Union{Element ,Vector{<:Element}}, b::Union{Element, Vector{<:Element}}; kwargs...
) = isapprox(transfer_matrix(a), transfer_matrix(b); kwargs...)

Base.isapprox(a::Matrix, b::Union{Element, Vector{<:Element}}; kwargs...) = isapprox(a, transfer_matrix(b); kwargs...)

Base.isapprox(a::Union{Element, Vector{<:Element}}, b::Matrix; kwargs...) = Base.isapprox(b, a) 



"""
    propagate(e::Union{Element, Vector{<:Element}}, b)

Propagate a beam `b` either by a single element `e::Element` or `Vector{<:Element}`.

Return is the final beam.
Also available as `e * b`.

## Example
```jldoctest
julia> beam = FreeSpace(10) * GeometricBeam(w=10.0, k=1.0)
GeometricBeam{Float64}(20.0, 1.0, 10.0)

julia> beam = [ThinLens(10), FreeSpace(10)] * GeometricBeam(w=10.0, k=0.0)
GeometricBeam{Float64}(0.0, -1.0, 10.0)

julia> beam.w ≈ 0
true
```
"""
function propagate(e::Element, b::GeometricBeam{T}) where T
    w, k = transfer_matrix(e) * [b.w, b.k]
    return GeometricBeam{T}(w, k, b.zpos + dz(e))
end

function propagate(e::Element, b::GaussianBeam{T}) where T
    ABCD = transfer_matrix(e)
    A, B, C, D = ABCD[1,1], ABCD[1,2], ABCD[2,1], ABCD[2,2]
    q_new = (A * b.q + B) / (C * b.q + D)
    return GaussianBeam{T}(q_new, b.zpos + dz(e), b.n, b.λ) 
end


function propagate(es::Vector{<:Element}, b::AbstractBeam)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

Base.:*(e::Union{Matrix, Element, Vector{<:Element}}, b::AbstractBeam) = propagate(e, b)
Base.:*(a::Element, b::Element) = transfer_matrix(a) * transfer_matrix(b)
Base.:*(a::Matrix, b::Element) = a * transfer_matrix(b)


"""
    trace(elems::Vector{<:Element}, b0::AbstractBeam)

Trace a beam `b0` through a vector of elements `elems`.
All intermediate states of the beam will be recorded.

Return is a `Vector` of states where the last entry is the final beam.
Final beam is equivalent to `propagate(elems, b0)`.


## Example
```jldoctest
julia> trace([ThinLens(10), FreeSpace(10)], GeometricBeam(w=10.0, k=0.0))
3-element Vector{GeometricBeam{Float64}}:
 GeometricBeam{Float64}(10.0, 0.0, 0.0)
 GeometricBeam{Float64}(10.0, -1.0, 0.0)
 GeometricBeam{Float64}(0.0, -1.0, 10.0)
```
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
