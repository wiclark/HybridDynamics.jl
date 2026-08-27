import HybridDynamics as HD
import Plots as plt

function f_circle(z, t)
    x, y, px, py = z
    return [px, py, 0, -1]
end
h_circle(z) = 1 - (z[1]^2 + z[2]^2)
function Δ_circle(z)
    x, y, px, py = z
    pxp = px - (2*x*(px*x + py*y))/(x^2 + y^2)
    pyp = py - (2*y*(px*x + py*y))/(x^2 + y^2)
    return [x, y, pxp, pyp]
end

sys = HD.GeneralSystem(f_circle, h_circle, Δ_circle; direction=-1)

function exponent_circle(x, y)
    if x^2+y^2 >= 1
        return NaN
    else
        return maximum(HD.LyapunovExponents(sys, [x, y, 0.0, 0.0]; run_length=10.0, run_iter=Int(1e3), transient=2000.0))
    end
end

res = 8
X = LinRange(-1, 1, res)
Y = LinRange(-1, 1, res)
Z = zeros(res, res)

@time Threads.@threads for i ∈ 1:res
    for j ∈ 1:res
        Z[i,j] = exponent_circle(X[i], Y[j])
    end
end

plt.heatmap(X, Y, Z', dpi=200)