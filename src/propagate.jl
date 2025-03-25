export propagate, trace

"""
    propagate(e::Union{Elements, Vector{<:Elements}}, b)

Propagate a beam `b` either by a single element `e::Element` or `e::Element3d` or `Vector{<:Elements}`.
If a Gaussian beams propagates through a 3D-element, it will become a 3D-Gaussian beam.

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
    return _n(T, e, b, q_new)
end

_n(::Type{T}, e::Element, b, q_new) where T = GaussianBeam{T}(q_new, b.zpos + dz(e), b.n, b.λ) 
_n(::Type{T}, e::Interface, b, q_new) where T = GaussianBeam{T}(q_new, b.zpos + dz(e), e.n2, b.λ) 

function propagate(e::Element3d, b::GaussianBeam{T}) where T
	ABCDsag, ABCDtan = transfer_matrix(e)
	Asag, Bsag, Csag, Dsag = ABCDsag[1,1], ABCDsag[1,2], ABCDsag[2,1], ABCDsag[2,2]
    qsag_new = (Asag * b.q + Bsag) / (Csag * b.q + Dsag)
    Atan, Btan, Ctan, Dtan = ABCDtan[1,1], ABCDtan[1,2], ABCDtan[2,1], ABCDtan[2,2]
    qtan_new = (Atan * b.q + Btan) / (Ctan * b.q + Dtan)
    return _n(T, e, b, qsag_new, qtan_new)
end

function propagate(e::Element3d, b::GaussianBeam3d{T}) where T
	ABCDsag, ABCDtan = transfer_matrix(e)
	Asag, Bsag, Csag, Dsag = ABCDsag[1,1], ABCDsag[1,2], ABCDsag[2,1], ABCDsag[2,2]
    qsag_new = (Asag * b.qsag + Bsag) / (Csag * b.qsag + Dsag)
    Atan, Btan, Ctan, Dtan = ABCDtan[1,1], ABCDtan[1,2], ABCDtan[2,1], ABCDtan[2,2]
    qtan_new = (Atan * b.qtan + Btan) / (Ctan * b.qtan + Dtan)
    return _n(T, e, b, qsag_new, qtan_new)
end

function propagate(e::Element, b::GaussianBeam3d{T}) where T
	ABCD = transfer_matrix(e)
	A, B, C, D = ABCD[1,1], ABCD[1,2], ABCD[2,1], ABCD[2,2]
    qsag_new = (A * b.qsag + B) / (C * b.qsag + D)
	qtan_new = (A * b.qtan + B) / (C * b.qtan + D)
    return _n(T, e, b, qsag_new, qtan_new)
end

_n(::Type{T}, e::Element3d, b::GaussianBeam, qsag_new, qtan_new) where T = GaussianBeam3d{T}(qsag_new, qtan_new, b.zpos + dz(e), b.n, b.λ)
_n(::Type{T}, e::Element3d, b::GaussianBeam3d, qsag_new, qtan_new) where T = GaussianBeam3d{T}(qsag_new, qtan_new, b.zpos + dz(e), b.n, b.λ)
_n(::Type{T}, e::Element, b::GaussianBeam3d, qsag_new, qtan_new) where T = GaussianBeam3d{T}(qsag_new, qtan_new, b.zpos + dz(e), b.n, b.λ)
_n(::Type{T}, e::Interface, b::GaussianBeam3d, q_new) where T = GaussianBeam3d{T}(qsag_new, qtan_new, b.zpos + dz(e), e.n2, b.λ)


function propagate(es::Vector{<:Elements}, b::AbstractBeam)
    return reduce((a,b) -> propagate(b, a), es, init=b)
end

Base.:*(e::Union{Matrix, Elements, Vector{<:Elements}}, b::AbstractBeam) = propagate(e, b)
Base.:*(a::Element, b::Element) = transfer_matrix(a) * transfer_matrix(b)
Base.:*(a::Elements, b::Elements) = transfer_matrix([a, b])
Base.:*(a::Matrix, b::Element) = a * transfer_matrix(b)
Base.:*(a::Matrix, b::Element3d) = a * transfer_matrix(b)
Base.:*(a::Matrix, b::Tuple{Matrix, Matrix}) = (a*b[1], a*b[2])
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

