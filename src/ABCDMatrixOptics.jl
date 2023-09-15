module ABCDMatrixOptics

using Parameters
using PrecompileTools


include("beam.jl")
include("elements.jl")
include("propagate.jl")
include("plots-recipes.jl")


@setup_workload begin
    system = [FreeSpace(100.0), FreeSpace(100), ThinLens(100), ThinLens(100.0), ThickLens(R1=10.0, R2=-1.0,t=1.0),
              Mirror(10), Mirror(10.9), Interface(n1=1.3, n2=1.4)] 
    @compile_workload begin
        system * GeometricBeam()
        system * GaussianBeam()
    end
end


end
