# Developer notes

Principles for every change to the cut-cell IBM code (2026-11-05).

1. **Numerical correctness is the highest priority.** It comes before agreement with a reference
   on a particular case, before speed and before convenience. A change has to be justified by the
   discretisation (consistency, order, conservation, independence of the decomposition), not only
   by a better number on one case.
2. **No case-specific fixes.** A fix that improves one case but does not follow from a general
   argument, or that has not been shown to hold across the validation suite, is not adopted. Agreement
   reached by cancelling errors is not an improvement; an option that helps one case and hurts
   another stays off by default (example: `ibmConvectionWallLaw`, `doc/secondOrderStudy.md`).
3. **Validation stresses one component at a time.** The validation cases are designed so that each
   component of the implementation (geometry, each operator, linking, wall treatment, wall law,
   turbulence model near walls, moving bodies, scalars, parallel decomposition) is tested as
   separately as possible. Complete separation is impossible; each case states which components it
   stresses and which it relies on, so that a failure points at a cause.
