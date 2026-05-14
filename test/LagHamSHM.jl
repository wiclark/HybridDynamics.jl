# testing

using HybridDynamics
using LaTeXStrings
import Plots as plt

init = [1.,0.];
tspan = (0., 4.);

function L(q,v)
	1/2 * v[1]^2 - q[1]^2
end;

prob = HybridDynamics.LagProb(L, init, tspan)
s = HybridDynamics.solve(prob, HybridDynamics.forward_euler; dt=0.01)

q = [x[1] for x in s[2]]
v = [x[2] for x in s[2]]
p1 = plt.plot(q, v, title = "Harmonic osc. phase portrait via Euler-Lagrange", label = "", xlabel = L"q", ylabel = L"\dot{q}")
display(p1)


function H(q,v)
	1/2 * v[1]^2 + q[1]^2
end;