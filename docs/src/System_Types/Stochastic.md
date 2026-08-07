### Stochastic

A stochastic system contains the data $(f, g, h, Δ; δ)$
```math
\begin{cases}
    dx = f(x,t)dt + g(x,t)dW, & h(x) \ne 0, \\
    x^+ = \Delta(x), & h(x) = 0
\end{cases}
```

```@autodocs
Modules = [HybridDynamics]
Pages = ["Systems/Stochastic.jl"]
```