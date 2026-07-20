# Changelog

## v0.4.0

### Testing

* **Filippov Testing:** Two new notebooks have been added as benchmarks for our Filippov system type. `FilippovSSF.ipynb` and `FilippovFuller.ipynb`. See those files for more details. 

* **ODE Solver Testing:** We have begun stress testing our integration methods. Look to `TestingNotebook.ipynb` for all our current results. Note Adaptive methods have not been fully tested yet. See `Known Bugs` below for more. 

### Known Bugs 
* **Adaptive Solvers:** `RK23` currently is bugged and does not work for most system types, especially `MechanicalSystems` and `FilippovSystems`. It can be used for simple examples but we recommend avoiding its use for more complex systems. It is also worth noting these integration methods (`RK23`, `RK45`, `AdaptiveABM2`, `AdaptiveABM3`) have not been fully tested rigorously. They will work but can be slightly off. Take their results with a grain of salt. 


## v0.3.0

### External Updates
* Added Documentation! See our dedicated website for the full list of functionality of this package.
* Updated logic of all solvers to utilize the current, predicted, and previous (from midpoint) states for use in quadratic event detection.
* All systems' `solve` dispatches now use `prob, solver; ...` instead of `prob; solver, ...` to better match `DifferentialEquations.jl`.
* Updated `pluto_examples` to reflect all our changes above. Give it a test!
* Linear and Affine systems now correctly identify if initial conditions start on the guard.
* Linear and Affine systems now have a toggle for logic order: Flow -> Event check, or Event check -> Flow. To utilize this, see the `event_before_flow` argument in the `solve` dispatch. 
* Made many different small tweaks to all solvers to allow (most of) them to work with all system types. See the Known Bugs section below for more details. 

### New Functions
* `extract_jumps`: Allows you to customize just the lines connecting discrete jumps within the plotting architecture. (See `split_jumps` to get rid of them completely!)

### Internal
* Unified variable names across all systems and solvers. This is mostly internal, but it may have changed some external calling (see Breaking Changes). 
* Added `take_step_systemtype!` to all system-specific `solve` dispatches. 

### Breaking Changes & Deprecations
* `jump_times` and `jump_indices` in Linear/Affine and General systems are now `event_times` and `event_indices` to match the rest of the package. 
* As `solve` dispatches are now `prob, solver; ...`, if you are calling `solve` in previous versions as `HybridDynamics.solve(prob, solver = ...)`, you will need to get rid of the `solver = ` keyword. Otherwise, Julia will silently default to the system type's default solver instead of using your specified one (usually `RK45` or `RK4`). 
* Let us know of any more we missed!

### Known Bugs & System Notes
* **General Systems:** `QuadraticLocator` and `BisectionLocator` (arguments for `EventLocator`) do not currently work for most solvers. They will run, but the output is not accurate. Currently, `RK45` is the only working solver for these.
* **Linear/Affine Systems:** Zeno detection is still primitive. Be wary of what it tells you and give us your findings so we can fine-tune! Be sure to utilize the tolerances in the `solve` dispatch to better fit your use case. 
* **Mechanical Systems:** Some results will 'fall out' of the accepted region. This is OK! It should only do this for a single step before catching itself and getting back on track. `RK23` struggles to provide accurate data. All event locators (except for Linear!) will not provide accurate data. Some are better than others, but be careful.
* **Nonholonomic Systems:** Some particularly bad solvers (such as Forward/Backwards Euler) will not behave nicely. `BDF2`, `ModifiedTrap`, `ModifiedMidpoint`, `RK23`, and `RK45` all provide inaccurate data with ALL locator options. 
* **Stochastic Systems:** All is in order (for now).
* **Filippov Systems:** Bug fixed to work with all adaptive solvers now!

---

## v0.2.0
* Finished implementing all system types. 

## v0.1.0
* Initial functionality.