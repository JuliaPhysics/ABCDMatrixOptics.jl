using ABCDMatrixOptics
using Test

@testset "ABCDMatrixOptics.jl" begin

    @testset "Base.≈" begin
        @test [ThinLens(100), ThinLens(200)] ≈ [ThinLens(200.0), ThinLens(100)]
    end

    @testset "Interface" begin
        b2 = (Interface(1.2) * GeometricBeam(x=1.0, k=1.0, n=1.0))
        @test b2 == GeometricBeam{Float64}(x=1.0, k=1.0 / 1.2, n=1.2, z = 0.0)
        @test Interface(1.2) * b2 == GeometricBeam{Float64}(x=1.0, k=1.0 / 1.2, n=1.2, z = 0.0)
        @test Interface(1.3) * b2 == GeometricBeam{Float64}(x=1.0, k=1.0 / 1.3, n=1.3, z = 0.0)
    end


    @testset "Beam trace" begin
        @test beamtrace([Interface(1.1), FreeSpace(100), ThinLens(100), FreeSpace(100)], GeometricBeam(x = 3.14)) == GeometricBeam{Float64}[GeometricBeam{Float64}(3.14, 0.0, 0.0, 1.0), GeometricBeam{Float64}(3.14, 0.0, 100.0, 1.0), GeometricBeam{Float64}(3.14, -0.031400000000000004, 100.0, 1.0), GeometricBeam{Float64}(-4.440892098500626e-16, -0.031400000000000004, 200.0, 1.0), GeometricBeam{Float64}(-4.440892098500626e-16, -0.028545454545454547, 200.0, 1.1)]


    end

    @testset "Free Space Propagation with Lens" begin
        f1 = FreeSpace(200)
        l1 = ThinLens(200.0)
        f12 = FreeSpace(200.0 + 300.0)
        l2 = ThinLens(300.0)
        f2 = FreeSpace(300.0)
        beam = GeometricBeam{Float64}(x=10.0, k=0.1)

        @test beam.x == 10.0
        @test beam.k == 0.1
        @test beam.n == 1.0
        @test beam.z == 0.0

        M = [f2, l2, f12, l1, f1]
        beam_p = propagate(M, beam)
        
        @test beam_p.x ≈ -15.0
        @test beam_p.k ≈ - 2/30
        @test beam_p.n == 1.0
        @test beam_p.z ≈ 1000.0

        beam_p2 = ABCDMatrixOptics.RTM(M) * [beam.x, beam.k]

        @test beam_p2 ≈ [-15.0, -0.06666666666666667] 
        @test beam_p2[1] ≈ beam_p.x
        @test beam_p2[2] ≈ beam_p.k
    end
end
