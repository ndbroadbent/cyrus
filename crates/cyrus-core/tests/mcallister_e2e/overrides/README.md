# Deprecated Overrides

These JSON files are legacy overrides from early reproduction work.
They are **not** used in first-principles runs.

Use the McAllister `.dat` files instead:

- `points.dat`, `heights.dat`, `dual_points.dat`, `dual_simplices.dat`
- `basis.dat`, `kklt_basis.dat`, `target_volumes.dat`
- `small_curves.dat`, `small_curves_gv.dat`

If you explicitly enable fixtures (`CYRUS_ALLOW_FIXTURES=1`), these files
may be used for comparison or debugging only. Do not rely on them for
production or validation.

**Note:** `CYRUS_ALLOW_FIXTURES` must not be used with `CYRUS_FIRST_PRINCIPLES`.
