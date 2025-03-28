module ABCDMatrixOpticsPlotExt


using ABCDMatrixOptics, Colors, Interpolations, RecipesBase
export beamplot3d

# mutable struct WithGaussianBeam3d
	# system::Vector{<:Elements}
	# beam::GaussianBeam3d
# end

mutable struct WithGaussianBeam
    system::Vector{<:Element}
    beam::GaussianBeam
end


mutable struct WithGeometricBeam
    system::Vector{<:Element}
    beam::GeometricBeam
end

# did not succed to plot 2 shape with a single @receipe ..
function Myplot3d(system::Vector{<:Elements}, beam::ABCDMatrixOptics.GaussianBeam3d)
    ds = ABCDMatrixOptics2.discretize(system, 200)
    N = length(ds) + 1
    Tz = typeof(float(beam.zpos))
    wsags = Vector{Tz}(undef, N)
    wtans = Vector{Tz}(undef, N)
    zs = Vector{Tz}(undef, N)
    wsags[1] = ABCDMatrixOptics.wsag(beam)
    wtans[1] = ABCDMatrixOptics.wtan(beam)
    zs[1] = beam.zpos
    beam = beam
    for i = 1:length(ds)
        beam = ds[i] * beam
        wsags[i+1] = ABCDMatrixOptics.wsag(beam)
        wtans[i+1] = ABCDMatrixOptics.wtan(beam)
        zs[i+1] = beam.zpos 
    end
    xs, ys_sag, ys_tan = vcat(zs, reverse(zs)), vcat(wsags, (-1.0) .* reverse(wsags)), vcat(wtans, (-1.0) .* reverse(wtans))
    shape_sag = [(xs[i], ys_sag[i], ) for i in 1:length(xs)]
    shape_tan = [(xs[i], ys_tan[i], ) for i in 1:length(xs)]
    index_iterfaces = findall(x -> !(x isa FreeSpace), ds)
    # We look for all non-FreeSpace elements to display them in the final plot
    z_interfaces = Vector{Tz}(undef, length(index_iterfaces))
    for i=1:length(z_interfaces)
        z_interfaces[i] = zs[index_iterfaces[i]]
    end
    plot(shape_sag, seriestype = :shape, 
            linecolor = color(beam.λ),
            fillcolor = color(beam.λ),
            fillalpha = 0.1,
            label = "Sag." * string(beam.λ),
            xlabel = "Distance along Beam Axis",
            yguide = "Extend along Beam Axis")
    plot!(shape_tan, seriestype = :shape,
            linecolor = color_tan(beam.λ),
            fillcolor = color_tan(beam.λ),
            fillalpha = 0.1,
            label = "Tan." * string(beam.λ))
    vline!(z_interfaces, primary = false, color = "black")
end

# type recipe, e.g. for `plot(WithBeam(system, beam))`
@recipe f(::Type{WithGaussianBeam}, data::WithGaussianBeam) = begin
    seriestype := :shape
    linecolor --> color(data.beam.λ)
    fillcolor --> color(data.beam.λ)
    fillalpha --> 0.1
    label --> string(data.beam.λ)
    xguide --> "Distance along Beam Axis"
    yguide --> "Extend along Beam Axis"
    # choose a very large number of discretizations to have some
    # chance of approximating minimum waist radii decently
    ds = ABCDMatrixOptics.discretize(data.system, 200)
    N = length(ds) + 1
    Tw = typeof(float(data.beam.zpos))
    Tz = typeof(float(data.beam.zpos))
    ws = Vector{Tw}(undef, N)
    zs = Vector{Tz}(undef, N)
    ws[1] = ABCDMatrixOptics.w(data.beam)
    zs[1] = data.beam.zpos
    beam = data.beam
    for i = 1:length(ds)
        beam = ds[i] * beam
        ws[i+1] = ABCDMatrixOptics.w(beam)
        zs[i+1] = beam.zpos 
    end
    xs, ys = vcat(zs, reverse(zs)), vcat(ws, (-1.0) .* reverse(ws))
    [(xs[i], ys[i]) for i in 1:length(xs)]
end


# type recipe, e.g. for `plot(WithBeam(system, beam))`
@recipe f(::Type{WithGeometricBeam}, data::WithGeometricBeam) = begin
    seriestype := :shape
    fillalpha --> 0.1
    label --> string("geometric beam")
    xguide --> "Distance along Beam Axis"
    yguide --> "Extend along Beam Axis"
    # choose a very large number of discretizations to have some
    # chance of approximating minimum waist radii decently
    ds = ABCDMatrixOptics.discretize(data.system, 200)
    N = length(ds) + 1
    Tw = typeof(float(data.beam.zpos))
    Tz = typeof(float(data.beam.zpos))
    ws = Vector{Tw}(undef, N)
    ws_mirror = Vector{Tw}(undef, N)
    zs = Vector{Tz}(undef, N)
    # for the geometric beam, trace beam wih negative k
    ws[1] = data.beam.w
    ws_mirror[1] = data.beam.w
    zs[1] = data.beam.zpos
    beam = data.beam
    beam_mirror = GeometricBeam(w=beam.w, zpos=beam.zpos, k=-beam.k)
    for i = 1:length(ds)
        beam = ds[i] * beam
        beam_mirror = ds[i] * beam_mirror
        ws[i+1] = beam.w
        ws_mirror[i+1] = beam_mirror.w
        zs[i+1] = beam.zpos 
    end
    xs, ys = vcat(zs, reverse(zs)), vcat(ws, reverse(ws_mirror))
    [(xs[i], ys[i]) for i in 1:length(xs)]
end




const _λ_hue_data = [
    380e-9 260 # "far" UV
    440e-9 250 # violet
    485e-9 180 # blue
    523e-9 120 # blue-green
    572e-9 60 # yellow
    620e-9 20
    700e-9 0 # (infra) red
    800e-9 0 # (infra) red
]
const hue_of_λ = Interpolations.LinearInterpolation(
    _λ_hue_data[:,1], _λ_hue_data[:,2]
)

"""

Return a (very approximately) correct color from the package
[Colors](http://docs.juliaplots.org/latest/colors/) for a given
wavelength (specified in SI units, i.e. in the range from
approximately `380e-9` to `700e-9`).

"""
function color(λ)
    # very rough transformation of λ to a color via package Colors
    λmin = 380e-9
    λmax = 800e-9
    if λ < λmin
        # ultraviolet (far spectral violet, dark)
        return Colors.HSL(260.0, 1.0, 0.25)
    elseif λ > λmax
        # infrared (red, dark)
        return Colors.HSL(0, 1.0, 0.25)
    end
    color = Colors.HSL(hue_of_λ(λ), 1.0, 0.5)
end

function color_tan(λ)
	#slightly darker color for tangential display to be different
    λmin = 380e-9
    λmax = 800e-9
    if λ < λmin
        return Colors.HSL(260.0, 1.0, 0.15)
    elseif λ > λmax
        return Colors.HSL(0, 1.0, 0.15)
    end
    color = Colors.HSL(hue_of_λ(λ), 1.0, 0.40)
end

# user recipe, e.g. for `plot(system, beam)`
# @recipe f(system::Vector{<:Elements}, beam::ABCDMatrixOptics.GaussianBeam3d) =
	# WithGaussianBeam3d(system, beam)
# @recipe f(system::Vector{<:Elements}, beam::ABCDMatrixOptics.GaussianBeam) =
	# WithGaussianBeam3d(system, beam)
@recipe f(system::Vector{<:Element}, beam::ABCDMatrixOptics.GaussianBeam) =
    WithGaussianBeam(system, beam)

@recipe f(system::Vector{<:Element}, beam::ABCDMatrixOptics.GeometricBeam) =
    WithGeometricBeam(system, beam)


end
