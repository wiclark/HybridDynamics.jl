### Mechanical

A mechanical system contains the data $(M, V, h)$ where $M$ is the (non-degenerate) mass matrix, $V$ is the potential engery, and $h$ is the event function. The problem is assumed to be in Hamiltonian form (the state is determined by position and momentum).
If a reset function is not specified, specular reflection will be assumed, i.e.,
```math
\Delta(p) = p - 2\frac{dh_q M(q)^{-1} p}{dh_q M(q)^{-1}dh_q} dh_q.
```

```@autodocs
Modules = [HybridDynamics]
Pages = ["Systems/Mechanical.jl"]
```