# Known Bugs

Here we document all known bugs and issues on the most current release.
---

* **Adaptive Solvers:** `RK23` currently is bugged and does not work for most system types, especially `MechanicalSystems` and `FilippovSystems`. It can be used for simple examples but we recommend avoiding its use for more complex systems. It is also worth noting these integration methods (`RK23`, `RK45`, `AdaptiveABM2`, `AdaptiveABM3`) have not been fully tested rigorously. They will work but can be slightly off. Take their results with a grain of salt. 