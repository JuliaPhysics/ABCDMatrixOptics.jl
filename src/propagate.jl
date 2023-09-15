export propagate, trace

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
Base.:*(a::Vector{<:Element}, b::Vector) = transfer_matrix(a) * b 


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
    bs = Vector{B}(undef, 1)
    bs[1] = b0
    for (idx, e) in enumerate(elems)
        push!(bs, trace(e, bs[end])...)
    end
    return bs
end


function trace(elems::Element, b0::B) where B
    return [propagate(elems, b0)]
end

