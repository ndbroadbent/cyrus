# Declared flux coordinate bases

`K_vec.dat` / `M_vec.dat` are expressed in the dual divisor basis chosen by
CYTools 2021 when the paper data was generated. That basis choice is part of
the declared model input (the flux integers are meaningless without it).
These files record it per example, extracted by running the vendored CYTools
2021 snapshot on each example's dual polytope + declared mirror chamber
(`string_theory/vendor/cytools_mcallister_2107`).

Pass to the runner via `--dual-basis <file>`. For 4-214-647 the computed
Cyrus dual basis happens to coincide with [3,4,5,8]'s transform behavior, so
the default works there; the other examples require the explicit declaration.
