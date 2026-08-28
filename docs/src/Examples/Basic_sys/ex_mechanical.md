## A Mechancal Example

!!! info "A Mechanical Problem"
	A Mechanical system contains the data $(M, V, G, N, R, E)$
	1. .$M(q)$ is the mass matrix.
	2. .$V(q)$ is the potential energy.
	3. .$G(q)$ is the event location function.
	4. .$N(q)$ is the normal to the guard.
	5. .$R(q)$ is the reset/impact map.
	6. .$E$ is the coefficient of restitution.
	The mechanical system has the Hamiltonian
	```math
		H(q,p) = \frac{1}{2}p^\top M^{-1}(q) p + V(q).
	```
	If the reset map is not specified, specular reflection is used.

```math
	H(x,y,p_x,p_y) = \frac{1}{2}(p_x^2+p_y^2) + y
```
```math
	h(x,y) = 1 - (x^2+y^2) \geq 0.
```

```@example mech
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings

M(q) = [1.0 0.0; 0.0 1.0]
V(q) = q[2]
h(q) = 1 - (q[1]^2+q[2]^2)
∇h(q) = [-2*q[1], -2*q[2]]
r = 0.6;
```

```@example mech
sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=r)
probM = HD.prob(sysM, [-0.95, 0.0, 0.2, -1.5], (0.0, 10.0))
solM = HD.solve(probM, HD.RK4())
```

```@example mech
xm = getindex.(solM.x, 1)
ym = getindex.(solM.x, 2)
plt.plot(xm, ym, label="Ball Trajectory", lw=2)
θ = LinRange(0, 2π, 100)
plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
```

The first time the state becomes Zeno can be found.

```@example mech
solM.zeno[1]
```