# Change Log

## v1.1.0

### New functionality
* New utility functions, `jump_interval` and `jump_count`, added. `jump_interval` takes an event index as input and returns a time interval where that event took place. `jump_count` has you input a time and outputs the amount of events that have happened at up to that time.
- New analysis function, 'LyapunovExponents()', added.
* New Solvers: `ExponentialSolver` has been added to the Linear and Affine system types. We have also added 3 new implicit solvers, `BackwardEuler`, `ImplicitTrap` and `RadauIIA`. 
* New "system" type. Added `GeneralTimeTriggered` which acts as a General system that inputs a set of times in place of the guard function. This allows you to indicate the times you want events to take place rather than utilizing a guard function. 

### Updates to old functionality
- Updated the Hermite interpolation's search for indicies around events and boundary conditions.

### Infrastructure
- Tests added for Filippov systems
- Various additions and improvements to documentation, including an "Analysis" section.


## v1.0.1

- Known bugs added to "Known Bugs" page.
- Various fixes to adaptive solvers. They now pass tests.
- Tests added for stochastic, and nonholonomic systems. All tests currently pass.
- Added check to see if the initial condition is on the guard to all system types. If so, the reset is applied immediately.
- Dense output via Hermite interpolation has been fixed for general, linear, affine, and mechanical systems, i.e. 'dx' is no longer off by one in its index. Attempting interpolation with a stochastic system now correctly throws an error.
- Many typos and small bugs (e.g. dropped negatives) have been fixed.


## v1.0.0

### Tests

* Tests have been added to the package. These will allow us to better keep up with how our updates work with the code as a whole. 
* Added `TestingNotebook.ipynb` which delves into detail on how we test our ODE solvers. Note for now it is quite sparse as we are focusing on actually making the solvers good but it will be updated in a future patch. This notebook includes many useful thigns so go take a look!

### ODE Solver Testing

* **Fixed Solvers:** Working much better now. They are not perfect currently but that will be fixed in the next update and they will work great for any applications of this package as of now. 
* **Adaptive Solvers:** Bugged at this time. Along with bringing them up to date with the new testing software (in `TestingNotebook.ipynb`) there has been some big behind the scenes changes for some of our solvers. `AdaptiveABM2` and `AdaptiveABM3` have been updated to use their adaptive coefficients as before we had fixed values. Then `RK23` and `RK45` have been given protection against catastrophic cancellation. Otherwise these solvers struggle at this time with certain systems, see the Known Bugs tab below. 

### Known Bugs

**Adaptive Solvers:** Adaptive solvers currently do not work as intended for `Mechanical`, `Filippov` and `General` systems. The fault lies in the more complicated guard structure. That in mind, the solvers should be more than reasonable for any non-hybrid applications as well as our other system types. These will be fixed in the next update. 

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