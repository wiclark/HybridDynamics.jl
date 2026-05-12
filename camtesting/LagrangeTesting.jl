# testing

using HybridDynamics
import Plots as plt

function L(q,v)
	1/2 * v[1]^2 - q[1]^2
end;

init = [1.,0.];
tspan = (0., 3.);

prob = HybridDynamics.LagProb(L, init, tspan)
s = HybridDynamics.solve(prob, HybridDynamics.forward_euler; dt=0.01)

q = [x[1] for x in s[2]]
v = [x[2] for x in s[2]]
plt.plot(q, v)