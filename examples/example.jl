### A Pluto.jl notebook ###
# v0.20.8

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
md"# Load the Package"

# ╔═╡ 82339a93-c98f-43e2-ba23-88c8c52bd45c
md"# Define Simple Elements"

# ╔═╡ e2e158c6-f690-4692-aae7-fecad5713b31
f1 = FreeSpace(200)

# ╔═╡ a4f123a0-7ff4-414d-b2ba-e997584a060f
l1 = ThinLens(200.0)

# ╔═╡ 7a19d973-7d69-45de-b731-2d673d89021b
l2 = ThinLens(100.0)

# ╔═╡ 02bd56b5-4fe1-4a70-a07d-1ffca857478c
f12 = FreeSpace(200 + 100)

# ╔═╡ b68ee752-9080-4bec-b0d0-9b65b6210aff
f2 = FreeSpace(100)

# ╔═╡ 5de0cb00-fa0b-4bb2-8f54-5abb4952c0eb
md"# Define a GeometricBeam"

# ╔═╡ a20eb02a-d12a-4815-87a2-2810f5f6cdeb
beam = GeometricBeam(w=10.0, k=0.1)

# ╔═╡ 8c916040-9d5e-4ec5-8e1b-b501faa2f911
beam.zpos

# ╔═╡ 566895eb-4471-4830-a2fc-3b3f77a18d54
M = [f2, l2, f12, l1, f1]

# ╔═╡ 1984aedc-e8a1-4331-a4ed-7ff3328ae348
md"# Propagate the beam through the system"

# ╔═╡ d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
M * beam

# ╔═╡ e67294e5-288a-4635-a139-3009b82e68f5
transfer_matrix(M) * [beam.w, beam.k]

# ╔═╡ 7ece0808-f64c-485f-bc23-91b244f474ab
md"# Store all intermediate states"

# ╔═╡ 0cb916dc-6057-46a6-8a37-945ee3e6d94d
trace(M, beam)

# ╔═╡ c174f8cb-a76d-4d55-8d07-88ef945ce4ab
md"# Plotting of Beams"

# ╔═╡ d4ba59e3-c8a9-4da2-ba7d-00a565a8fc80
Pkg.instantiate()

# ╔═╡ c1dbc6a7-f302-45d0-bcc8-9e5044bc254b
red_beam = GaussianBeam(w0=5e-3)

# ╔═╡ 8eb11a1c-65d1-4d5a-aeae-d9d1e6c6a988
blue_beam = GaussianBeam(w0=5e-3, λ=405e-9)

# ╔═╡ 36184904-8683-4f9f-968e-b56ffe7d0428
beam_geometric = GeometricBeam(w=5e-3, k=0.1e-3)

# ╔═╡ 8faab7e3-56a9-4218-b187-7039dd43cb7b
md"
We see that the 4f system both magnifies the red and blue beam by a factor of 2.

But, the intermediate beam size and the curvatures depend on the wavelength!
"

# ╔═╡ 5844cbf9-1551-4124-a168-5c7e1d4e4746
Revise.retry()

# ╔═╡ 332416d5-91b8-497b-8d5e-57d8e1bc9a4d
begin
	plot(M, red_beam)
	plot!(M, blue_beam)
	plot!(M, beam_geometric)
end

# ╔═╡ 720bc546-8196-4210-9cb2-364799951ccf
md"# Plot Simple Lens system"

# ╔═╡ 1a9aef8d-80a0-4dc2-9f63-3bd9ee52d58d
tl1 = ThickLens(R1=100, R2=-100, t=0)

# ╔═╡ 3832f547-0d24-4688-b77d-453053e577ea
tl2 = ThickLens(R1=200, R2=-200, t=50)

# ╔═╡ 38ecda12-6e81-4035-ad40-9498334471e7
tl1.focal_length

# ╔═╡ 8579099b-9c7e-47b8-8114-5c90653e2df0
tl2.focal_length

# ╔═╡ 3442fbeb-99a4-4af7-9816-283b533c3421
M2 = [f2, tl1, f12, tl2, f1]

# ╔═╡ ad8215d3-901e-4629-98c0-ee56533e5c1c
begin
	plot(M2, beam)
	# specify the height and it'll show the lens with that height
	plot!(M2, height=40)
end

# ╔═╡ 51c61aa0-4e1c-4c22-8c4f-8b6611ed5075
md"# Effects of refractive index on Beam"

# ╔═╡ 829c357a-f490-4ad8-98ce-b4ebbf09c552
plot([FreeSpace(10e-3), Interface(n1=1.0, n2=2.0), FreeSpace(10e-3)], GaussianBeam(w0=10f-6))

# ╔═╡ fc8c7aba-6a28-462f-a634-e25ead9486c0


# ╔═╡ d6b93543-d91b-490a-8b21-2f04f7b781d3
md"# ThickLensTracing"

# ╔═╡ c1a34ab4-f75d-4dda-bd38-d3f0c22ee793
system2 = [FreeSpace(10e-3), ThickLens(R1=10e-3, R2=-10e-3, t=10e-3), FreeSpace(10e-3)]

# ╔═╡ d1dec156-89ea-477f-88c5-069bc73e9e6b
b2 = GaussianBeam(w0=40e-6)

# ╔═╡ 6da1508c-5bd9-4278-bde0-d894d2b932f2
trace(system2, b2)

# ╔═╡ 94a7cb0e-1ea4-42d7-833b-2a92342451dc
plot(system2, b2)

# ╔═╡ Cell order:
# ╠═73c9f846-50de-11ee-1b1b-2f253687c2e6
# ╟─c91dcb4a-73f6-42b3-a930-a21965c8c9dc
# ╟─6cd52e19-8721-45a4-8fb1-a96b4ad004f1
# ╠═1a4560a3-da99-4ee5-97f9-ab2d97e641da
# ╟─82339a93-c98f-43e2-ba23-88c8c52bd45c
# ╠═e2e158c6-f690-4692-aae7-fecad5713b31
# ╠═a4f123a0-7ff4-414d-b2ba-e997584a060f
# ╠═7a19d973-7d69-45de-b731-2d673d89021b
# ╠═02bd56b5-4fe1-4a70-a07d-1ffca857478c
# ╠═b68ee752-9080-4bec-b0d0-9b65b6210aff
# ╟─5de0cb00-fa0b-4bb2-8f54-5abb4952c0eb
# ╠═a20eb02a-d12a-4815-87a2-2810f5f6cdeb
# ╠═8c916040-9d5e-4ec5-8e1b-b501faa2f911
# ╠═566895eb-4471-4830-a2fc-3b3f77a18d54
# ╟─1984aedc-e8a1-4331-a4ed-7ff3328ae348
# ╠═d58bcfb6-6c7c-4870-9d30-c0093bc6dea8
# ╠═e67294e5-288a-4635-a139-3009b82e68f5
# ╟─7ece0808-f64c-485f-bc23-91b244f474ab
# ╠═0cb916dc-6057-46a6-8a37-945ee3e6d94d
# ╟─c174f8cb-a76d-4d55-8d07-88ef945ce4ab
# ╠═6389232e-ed13-40e1-8122-9875156c0603
# ╠═d4ba59e3-c8a9-4da2-ba7d-00a565a8fc80
# ╠═c1dbc6a7-f302-45d0-bcc8-9e5044bc254b
# ╠═8eb11a1c-65d1-4d5a-aeae-d9d1e6c6a988
# ╠═36184904-8683-4f9f-968e-b56ffe7d0428
# ╟─8faab7e3-56a9-4218-b187-7039dd43cb7b
# ╠═5844cbf9-1551-4124-a168-5c7e1d4e4746
# ╠═332416d5-91b8-497b-8d5e-57d8e1bc9a4d
# ╟─720bc546-8196-4210-9cb2-364799951ccf
# ╠═1a9aef8d-80a0-4dc2-9f63-3bd9ee52d58d
# ╠═3832f547-0d24-4688-b77d-453053e577ea
# ╠═38ecda12-6e81-4035-ad40-9498334471e7
# ╠═8579099b-9c7e-47b8-8114-5c90653e2df0
# ╠═3442fbeb-99a4-4af7-9816-283b533c3421
# ╠═ad8215d3-901e-4629-98c0-ee56533e5c1c
# ╟─51c61aa0-4e1c-4c22-8c4f-8b6611ed5075
# ╠═829c357a-f490-4ad8-98ce-b4ebbf09c552
# ╠═fc8c7aba-6a28-462f-a634-e25ead9486c0
# ╟─d6b93543-d91b-490a-8b21-2f04f7b781d3
# ╠═c1a34ab4-f75d-4dda-bd38-d3f0c22ee793
# ╠═d1dec156-89ea-477f-88c5-069bc73e9e6b
# ╠═6da1508c-5bd9-4278-bde0-d894d2b932f2
# ╠═94a7cb0e-1ea4-42d7-833b-2a92342451dc
