### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 6dab3570-7290-11f1-b67c-277aa2a67d52
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ 2b7869d2-dcd8-4ccf-b088-605b0c261b0b
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 78e500a7-75da-44f3-adfc-5c649de68608
function f_ball(x, t)
	g = 9.81
	α = 0.0
	q, v = x
	return [v, -g-α*v*abs(v)]
end

# ╔═╡ f3bcbae6-5f47-45d7-be2b-d44c166101bb
h_ball(x) = x[1]

# ╔═╡ 838e3bc0-c0b5-416f-9cef-ba1f1d44e450
Δ_ball(x) = [x[1], -0.8*x[2]]

# ╔═╡ d329682f-fb7b-49d6-8bed-b45f83109924
Δ_ball([8.600118115297158e-7, 0.718907768224883])

# ╔═╡ b8befa98-6049-424b-8f4c-e7d868dea35d
begin
	sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
	probG = HD.prob(sysG, [10.0, 0.0], (0.0, 15.0))
	solG = HD.solve(probG)
end

# ╔═╡ 262d6990-173d-4f11-92bd-647874bc4c11
solG.x

# ╔═╡ 338eeb44-86d5-4827-bb60-488b12269924
function get_rid_of_jump_lines(sol)
    t = Float64[]
    x1 = Float64[]
    x2 = Float64[]

    for i in 1:length(sol.x)
        push!(t, sol.t[i])
        push!(x1, sol.x[i][1])
        push!(x2, sol.x[i][2])

        # If next timestamp is identical (jump event), break the line
        if i < length(sol.t) && sol.t[i] == sol.t[i+1]
            push!(t, NaN)
            push!(x1, NaN)
            push!(x2, NaN)
        end
    end
    return t, x1, x2
end

# ╔═╡ 49da0e79-d4ec-4dc4-b14c-206f08f9f0bc
begin
	t_plot, x1_plot, x2_plot = get_rid_of_jump_lines(solG)
	plt.plot(t_plot, x1_plot, label="Position", lw=2, lc=:blue)
	plt.plot!(t_plot, x2_plot, label="Velocity", lw=2, color=:orange)
	plt.plot!(title = "Nonlinear bouncing ball", xlabel = "Time", grid = true)
end

# ╔═╡ 13da5896-c5ef-4829-b934-917a42495d01
begin
	M(q) = [1.0]
	V(q) = 9.81*q[1]
	h(q) = q[1]
	∇h(q) = [1.0]
end

# ╔═╡ 88bd4e41-04e8-4fa2-9223-66b25b5bd6b0
begin
	sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=0.8)
	probM = HD.prob(sysM, [10.0, 0.0], (0.0, 15.0))
	solM = HD.solve(probM, solver=HD.RK4())
end

# ╔═╡ c644ac64-4bc5-45ca-9aa2-60f0170e5eac
plt.plot(solM.t, getindex.(solM.x, 1))

# ╔═╡ 16fad695-ebd3-4470-a91a-96d41ae202f6
begin
	plt.plot(t_plot, x1_plot, label="Position", lw=3, lc=:blue)
	plt.plot!(solM.t, getindex.(solM.x, 1))
end

# ╔═╡ Cell order:
# ╠═6dab3570-7290-11f1-b67c-277aa2a67d52
# ╠═2b7869d2-dcd8-4ccf-b088-605b0c261b0b
# ╠═78e500a7-75da-44f3-adfc-5c649de68608
# ╠═f3bcbae6-5f47-45d7-be2b-d44c166101bb
# ╠═838e3bc0-c0b5-416f-9cef-ba1f1d44e450
# ╠═d329682f-fb7b-49d6-8bed-b45f83109924
# ╠═b8befa98-6049-424b-8f4c-e7d868dea35d
# ╠═49da0e79-d4ec-4dc4-b14c-206f08f9f0bc
# ╠═262d6990-173d-4f11-92bd-647874bc4c11
# ╟─338eeb44-86d5-4827-bb60-488b12269924
# ╠═13da5896-c5ef-4829-b934-917a42495d01
# ╠═88bd4e41-04e8-4fa2-9223-66b25b5bd6b0
# ╠═c644ac64-4bc5-45ca-9aa2-60f0170e5eac
# ╠═16fad695-ebd3-4470-a91a-96d41ae202f6
