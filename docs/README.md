# Cyrus Docs

Current local notes:

- [CYGV_AUDIT.md](CYGV_AUDIT.md): CYTools and `cygv` source audit for the
  McAllister GV pipeline, including the split between low-dimensional generic
  `compute_gvs()` and high-dimensional selected toric curves.
- [CKYZ_SERIES_DOMAIN_AUDIT.md](CKYZ_SERIES_DOMAIN_AUDIT.md): focused source
  audit for the finite monomial/coefficient domain needed by local CKYZ
  potent-ray rows.
- [GA_READY_PIPELINE_AUDIT.md](GA_READY_PIPELINE_AUDIT.md): checklist mapping
  the GA-ready objective to current code/test evidence and unresolved blockers.
- [MCALLISTER_DATA_POLICY.md](MCALLISTER_DATA_POLICY.md): allowed ancillary
  inputs, validation checkpoints, replay-only artifacts, and current unresolved
  gaps for the no-replay McAllister runner.
- [POTENT_RAY_SOURCE_READ.md](POTENT_RAY_SOURCE_READ.md): potent-ray source
  audit, current rank-two/reflexive-polygon inventory, and no-cheat boundary
  for `potent_rays*.dat`.
- [LOCAL_TORIC_GV_SOURCE_MAP.md](LOCAL_TORIC_GV_SOURCE_MAP.md): source map for
  the next non-`P^2` local toric GV implementation, centered on CKYZ local
  mirror symmetry and the topological vertex fallback.
- [RUST_CRATES_GUIDE.md](RUST_CRATES_GUIDE.md): notes on Rust dependency usage.

For broader project-level formula and architecture notes, see
`project/project_docs/`.
