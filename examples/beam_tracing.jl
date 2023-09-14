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
using ABCDMatrixOptics, PlutoUI

# ╔═╡ 6389232e-ed13-40e1-8122-9875156c0603
using Plots

# ╔═╡ c91dcb4a-73f6-42b3-a930-a21965c8c9dc
TableOfContents()

# ╔═╡ 6cd52e19-8721-45a4-8fb1-a96b4ad004f1
md"## Load the Package"

# ╔═╡ 82339a93-c98f-43e2-ba23-88c8c52bd45c
md"## Define Simple Elements"

# ╔═╡ e2e158c6-f690-4692-aae7-fecad5713b31
f1 = FreeSpace(200)

# ╔═╡ a4f123a0-7ff4-414d-b2ba-e997584a060f
l1 = ThinLens(200.0)

# ╔═╡ 02bd56b5-4fe1-4a70-a07d-1ffca857478c
f12 = FreeSpace(200 + 100)

# ╔═╡ 38603793-c1f2-4281-bb7d-874ccdab6fe0
l2 = ThinLens(100.0)

# ╔═╡ b68ee752-9080-4bec-b0d0-9b65b6210aff
f2 = FreeSpace(100)

# ╔═╡ 5de0cb00-fa0b-4bb2-8f54-5abb4952c0eb
md"## Define a GeometricBeam"

# ╔═╡ a20eb02a-d12a-4815-87a2-2810f5f6cdeb
beam = GeometricBeam{Float64}(w=10.0, k=0.1)

# ╔═╡ 566895eb-4471-4830-a2fc-3b3f77a18d54
M = [f2, l2, f12, l1, f1]

# ╔═╡ 1984aedc-e8a1-4331-a4ed-7ff3328ae348
md"## Propagate the beam through the system"

# ╔═╡ d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
M * beam

# ╔═╡ e67294e5-288a-4635-a139-3009b82e68f5
transfer_matrix(M) * [beam.w, beam.k]

# ╔═╡ 7ece0808-f64c-485f-bc23-91b244f474ab
md"## Store all intermediate states"

# ╔═╡ 0cb916dc-6057-46a6-8a37-945ee3e6d94d
trace(M, beam)

# ╔═╡ c174f8cb-a76d-4d55-8d07-88ef945ce4ab
md"## Plotting of Gaussian Beams"

# ╔═╡ c1dbc6a7-f302-45d0-bcc8-9e5044bc254b
red_beam = GaussianBeam(w0=5e-3)

# ╔═╡ 8eb11a1c-65d1-4d5a-aeae-d9d1e6c6a988
blue_beam = GaussianBeam(w0=5e-3, λ=405e-9)

# ╔═╡ 8faab7e3-56a9-4218-b187-7039dd43cb7b
md"
We see that the 4f system both magnifies the red and blue beam by a factor of 2.

But, the intermediate beam size and the curvatures depend on the wavelength!
"

# ╔═╡ 332416d5-91b8-497b-8d5e-57d8e1bc9a4d
begin
	plot(M, red_beam)
	plot!(M, blue_beam)
end

# ╔═╡ Cell order:
# ╠═73c9f846-50de-11ee-1b1b-2f253687c2e6
# ╟─c91dcb4a-73f6-42b3-a930-a21965c8c9dc
# ╟─6cd52e19-8721-45a4-8fb1-a96b4ad004f1
# ╠═1a4560a3-da99-4ee5-97f9-ab2d97e641da
# ╟─82339a93-c98f-43e2-ba23-88c8c52bd45c
# ╠═e2e158c6-f690-4692-aae7-fecad5713b31
# ╠═a4f123a0-7ff4-414d-b2ba-e997584a060f
# ╠═02bd56b5-4fe1-4a70-a07d-1ffca857478c
# ╠═38603793-c1f2-4281-bb7d-874ccdab6fe0
# ╠═b68ee752-9080-4bec-b0d0-9b65b6210aff
# ╟─5de0cb00-fa0b-4bb2-8f54-5abb4952c0eb
# ╠═a20eb02a-d12a-4815-87a2-2810f5f6cdeb
# ╠═566895eb-4471-4830-a2fc-3b3f77a18d54
# ╟─1984aedc-e8a1-4331-a4ed-7ff3328ae348
# ╠═d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
# ╠═e67294e5-288a-4635-a139-3009b82e68f5
# ╟─7ece0808-f64c-485f-bc23-91b244f474ab
# ╠═0cb916dc-6057-46a6-8a37-945ee3e6d94d
# ╟─c174f8cb-a76d-4d55-8d07-88ef945ce4ab
# ╠═6389232e-ed13-40e1-8122-9875156c0603
# ╠═c1dbc6a7-f302-45d0-bcc8-9e5044bc254b
# ╠═8eb11a1c-65d1-4d5a-aeae-d9d1e6c6a988
# ╟─8faab7e3-56a9-4218-b187-7039dd43cb7b
# ╠═332416d5-91b8-497b-8d5e-57d8e1bc9a4d
