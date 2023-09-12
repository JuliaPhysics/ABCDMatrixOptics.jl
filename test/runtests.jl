using ABCDMatrixOptics
using Test

@testset "ABCDMatrixOptics.jl" begin
    f1 = FreeSpace(200)
    l1 = ThinLens(200.0)
    f12 = FreeSpace(200 + 300)
    l2 = ThinLens(300.0)
    f2 = FreeSpace(300)
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

    beam_p2 = RTM(M) * [beam.x, beam.k]

    @test beam_p2 ≈ [-15.0, -0.06666666666666667] 
    @test beam_p2[1] ≈ beam_p.x
    @test beam_p2[2] ≈ beam_p.k
end
