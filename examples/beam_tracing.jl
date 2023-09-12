### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 73c9f846-50de-11ee-1b1b-2f253687c2e6
begin
	using Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	using Revise
end

# ╔═╡ 1a4560a3-da99-4ee5-97f9-ab2d97e641da
using ABCDMatrixOptics#:w, Plots

# ╔═╡ e2e158c6-f690-4692-aae7-fecad5713b31
f1 = FreeSpace(200)

# ╔═╡ a4f123a0-7ff4-414d-b2ba-e997584a060f
l1 = ThinLens(200.0)

# ╔═╡ 02bd56b5-4fe1-4a70-a07d-1ffca857478c
f12 = FreeSpace(200 + 300)

# ╔═╡ 38603793-c1f2-4281-bb7d-874ccdab6fe0
l2 = ThinLens(300.0)

# ╔═╡ b68ee752-9080-4bec-b0d0-9b65b6210aff
f2 = FreeSpace(300)

# ╔═╡ a20eb02a-d12a-4815-87a2-2810f5f6cdeb
@show beam = GeometricBeam{Float64}(x=10.0, k=0.1)

# ╔═╡ 566895eb-4471-4830-a2fc-3b3f77a18d54
M = [f2, l2, f12, l1, f1]

# ╔═╡ aac37a72-ca16-4c89-ab26-e484ca5b3632
Revise.retry()

# ╔═╡ d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
propagate(M, beam)

# ╔═╡ d3935d82-1cda-4b47-8baa-df244fba389f
propagate(M, beam)

# ╔═╡ e67294e5-288a-4635-a139-3009b82e68f5
@show RTM(M) * [beam.x, beam.k]

# ╔═╡ 0cb916dc-6057-46a6-8a37-945ee3e6d94d
beamtrace(M, beam)

# ╔═╡ Cell order:
# ╠═73c9f846-50de-11ee-1b1b-2f253687c2e6
# ╠═1a4560a3-da99-4ee5-97f9-ab2d97e641da
# ╠═e2e158c6-f690-4692-aae7-fecad5713b31
# ╠═a4f123a0-7ff4-414d-b2ba-e997584a060f
# ╠═02bd56b5-4fe1-4a70-a07d-1ffca857478c
# ╠═38603793-c1f2-4281-bb7d-874ccdab6fe0
# ╠═b68ee752-9080-4bec-b0d0-9b65b6210aff
# ╠═a20eb02a-d12a-4815-87a2-2810f5f6cdeb
# ╠═566895eb-4471-4830-a2fc-3b3f77a18d54
# ╠═aac37a72-ca16-4c89-ab26-e484ca5b3632
# ╠═d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
# ╠═d3935d82-1cda-4b47-8baa-df244fba389f
# ╠═e67294e5-288a-4635-a139-3009b82e68f5
# ╠═0cb916dc-6057-46a6-8a37-945ee3e6d94d
