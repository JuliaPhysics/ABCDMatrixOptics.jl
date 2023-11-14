using ABCDMatrixOptics
using Test
using Plots

@testset "ABCDMatrixOptics.jl" begin

    @testset "Base.≈" begin
        @test [ThinLens(100), ThinLens(200)] ≈ [ThinLens(200.0), ThinLens(100)]
    end

    @testset "Interface" begin
        b2 = (Interface(n1=1.0, n2=1.2) * GeometricBeam(w=1.0, k=1.0))
        @test b2 == GeometricBeam(w=1.0, k=1.0 / 1.2, zpos = 0.0)
        @test Interface(n1=1.2, n2=1.2) * b2 == GeometricBeam{Float64}(1.0, 1.0 / 1.2, 0.0)
        @test Interface(n1=1.2, n2=1.3) * b2 == GeometricBeam{Float64}(1.0, 1.0 / 1.3, 0.0)

        @test ABCDMatrixOptics.dz([1 0; 0 1]) == Inf

        @test Interface(1.0, 1) == Interface{Float64}(1.0, 1.0, Inf)
        @test Interface(1, 2) == Interface{Float64}(1.0, 2.0, Inf)
        @test Interface(1, 2, 1.0) == Interface{Float64}(1.0, 2.0, 1.0)
        @test Interface(1.0f0, 2.0f0, 1.0f0) == Interface{Float32}(1.0f0, 2.0f0, 1.0f0)
    end

    @testset "ThickLens" begin
        @test ≈(([ThickLens(R1=100.0, R2=-50.0, t=20.0, n_lens=1.3), FreeSpace(110)]* GeometricBeam(w=1.0, k=0.0)).w + 1, 1, rtol=0.007)
        @test transfer_matrix(ThickLens(R1=-12.1, R2=20.0, t=31.1, n_lens=1.3, n1=2.0, n2=1.4)) ≈
        [1 0; (1.3 - 1.4)/ 20 / 1.4 1.3 / 1.4] * [1 31.1; 0 1] *  [1 0; (2 - 1.3) / 1.3 / (-12.1) 2 / 1.3]

        tr = trace([ThickLens(R1=100.0, R2=-50.0, t=20.0, n_lens=1.3), FreeSpace(110)], GeometricBeam(w=1.0, k=0.0))

        @test tr == GeometricBeam{Float64}[GeometricBeam{Float64}(1.0, 0.0, 0.0), GeometricBeam{Float64}(0.9538461538461538, -0.008723076923076924, 20.0), GeometricBeam{Float64}(-0.005692307692307774, -0.008723076923076924, 130.0)]
        @test ThickLens(R1=-12.1, R2=20.0, t=0.0, n_lens=1.3) ≈ ThinLens(f=inv((1.3 - 1) * (-1 / 12.1 - 1 / 20))) 
    end


    @testset "Mirror" begin
        @test [FreeSpace(100), Mirror(-100), FreeSpace(-50)] * GeometricBeam(w=100.0, k=0.0) == GeometricBeam{Float64}(0.0, 2.0, 50.0)

    end


    @testset "Propagate" begin
        b0 = [Interface(n1=1.0, n2=1.2), FreeSpace(100)] * GeometricBeam(w = 1.0, k = 1.0)
        @test b0 == GeometricBeam{Float64}(84.33333333333334, 0.8333333333333334, 100.0)

        b1 = ABCDMatrixOptics.transfer_matrix([Interface(n1=1.0, n2=1.2), FreeSpace(100)]) * [1.0, 1.0]

        @test ABCDMatrixOptics.transfer_matrix([1 0; -1 -1]) == [1 0; -1 -1]
        @test [b0.w, b0.k] == b1
        @test [ThinLens(100), FreeSpace(100)] * [100, 0.0] ≈ [0.0, -1.0]
    end

    @testset "ThinLens" begin
        @test ThinLens(-100, 100, 3, 1.3).f ≈ inv((3-1.3) / 1.3 * (1/(-100) - 1/100))
    end


    @testset "Beam trace" begin
        @test trace(reverse([Interface(n1=1.0, n2=1.1), FreeSpace(100), ThinLens(100), FreeSpace(100)]), GeometricBeam(w = 3.14)) == GeometricBeam{Float64}[GeometricBeam{Float64}(3.14, 0.0, 0.0), GeometricBeam{Float64}(3.14, 0.0, 100.0), GeometricBeam{Float64}(3.14, -0.031400000000000004, 100.0), GeometricBeam{Float64}(-4.440892098500626e-16, -0.031400000000000004, 200.0), GeometricBeam{Float64}(-4.440892098500626e-16, -0.028545454545454547, 200.0)]


    end

    @testset "Free Space Propagation with Lens" begin
        f1 = FreeSpace(200)
        l1 = ThinLens(200.0)
        f12 = FreeSpace(200.0 + 300.0)
        l2 = ThinLens(300.0)
        f2 = FreeSpace(300.0)
        beam = GeometricBeam(w=10.0, k=0.1)

        @test beam.w == 10.0
        @test beam.k == 0.1
        @test beam.zpos == 0.0

        M = [f1, l1, f12, l2, f2]
        beam_p = propagate(M, beam)
        
        @test beam_p.w ≈ -15.0
        @test beam_p.k ≈ - 2/30
        @test beam_p.zpos ≈ 1000.0

        beam_p2 = ABCDMatrixOptics.transfer_matrix(M) * [beam.w, beam.k]

        @test beam_p2 ≈ [-15.0, -0.06666666666666667] 
        @test beam_p2[1] ≈ beam_p.w
        @test beam_p2[2] ≈ beam_p.k
        @test f2 * l2 * f12 * l1 * f1 ≈ M
        @test M ≈ f2 * l2 * f12 * l1 * f1
    end

    @testset "discretize" begin
        @test ABCDMatrixOptics.discretize(FreeSpace(100), 10) == fill(FreeSpace(10.0), 10)
        @test ABCDMatrixOptics.discretize(ThinLens(100), 10) == ThinLens(100)
        @test ABCDMatrixOptics.discretize([FreeSpace(100), FreeSpace(2)], 2) == [FreeSpace(50.0), FreeSpace(50.0), FreeSpace(1.0), FreeSpace(1.0)]


    end

    @testset "Gaussian Beam" begin
        beam = GaussianBeam(w0=100e-6)
        @test ABCDMatrixOptics.zR(beam) ≈ 0.049630215696521214
        @test beam.q ≈ 0.0 + 0.049630215696521214im
        
        beam2 = FreeSpace(100e-3) * GaussianBeam(w0=100e-6)
        @test beam2.q ≈ 0.1 + 0.049630215696521214im
        @test ABCDMatrixOptics.R(beam2) ≈ 0.12463158310083221
        @test ABCDMatrixOptics.zR(beam2) ≈ 0.049630215696521214
        @test ABCDMatrixOptics.w(beam2) ≈ 0.00022494062272623123


        beam3 = FreeSpace(100e-3) * GaussianBeam(w0=100e-6, λ=100e-9, n=1.3, zpos=0)
        @test beam3.q ≈ GaussianBeam(100e-3 + 1im * π * 1.3 * 100e-6^2 / 100e-9).q 

        beam4 = Interface(n1=1.0, n2=2.0) * GaussianBeam()
        @test beam4.n == 2.0
        @test ABCDMatrixOptics.zR(beam4) /2 ≈ ABCDMatrixOptics.zR(GaussianBeam())
    end

    # how to do that properly?
    @testset "Plots" begin
        p = plot([FreeSpace(100e-3)], GaussianBeam(w0=100e-6, λ=500e-9, n=1.3, zpos=0))
        @test p == p
        p = plot([FreeSpace(100e-3)], GaussianBeam(w0=100e-6, λ=1000e-9, n=1.3, zpos=0))
        @test p == p
        p = plot([FreeSpace(100e-3)], GaussianBeam(w0=100e-6, λ=100e-9, n=1.3, zpos=0))
        @test p == p
    end

    return true
end
