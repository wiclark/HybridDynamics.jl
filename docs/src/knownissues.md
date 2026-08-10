# Known Bugs

Here we document all known bugs and issues on the most current release.
---

## v1.0.1 Known Bugs

* **General Systems:** Currently the `GeneralSystem` type struggles with a few of our solvers. Testing has shown the Fixed LMM solvers (`AdamsBashforth2`, `AdamsBashforth3` and `BDF2`) all fall short of our expected results when utilized with General systems. It is also known that `RK45()` struggles to work with General Systems exceptional pathology catcher, along with showing the same issues as the LMM solvers above. 
* **RK45:** While `RK45()` now works unlike the last version, it still struggles with any Zeno-like behavior.
* **Hermite interpolation:** Nonholonomic and Filippov systems is incorrect when using Hermite Interpolation.

