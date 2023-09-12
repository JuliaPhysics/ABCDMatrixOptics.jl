mutable struct WithBeam
    system::Vector{<:Element}
    beam::AbstractBeam
end

# type recipe, e.g. for `plot(WithBeam(system, beam))`
@recipe f(::Type{WithBeam}, data::WithBeam) = begin
    seriestype := :shape
 #linecolor --> color(data.beam.λ)
 #   fillcolor --> color(data.beam.λ)
    fillalpha --> 0.1
 #  label --> string(data.beam.λ)
    aspect_ratio --> :equal
    xlabel --> "Distance along Beam Axis"
    ylabel --> "1/e^2 Boundary"
    # choose a very large number of discretizations to have some
    # chance of approximating minimum waist radii decently
    ds = discretize(data.system, 200)
    N = length(ds) + 1
    Tw = typeof(float(spotsize(data.beam)))
    Tz = typeof(float(location(data.beam)))
    ws = Vector{Tw}(undef, N)
    zs = Vector{Tz}(undef, N)
    ws[1] = spotsize(data.beam)
    zs[1] = location(data.beam)
    beam = data.beam
    for i = 1:length(ds)
        beam = transform(ds[i], beam)
        ws[i+1] = spotsize(beam)
        zs[i+1] = location(beam)
    end
    xs, ys = vcat(zs, reverse(zs)), vcat(ws, (-1.0) .* reverse(ws))
    [(xs[i], ys[i]) for i in 1:length(xs)]
end

# user recipe, e.g. for `plot(system, beam)`
@recipe f(system::Vector{<:Element}, beam::AbstractBeam) =
    WithBeam(system, beam)
