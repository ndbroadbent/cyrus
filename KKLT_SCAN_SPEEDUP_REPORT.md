# Speeding up KKLT vacuum verification: a basis-selection / chamber-search bottleneck

**Audience:** external research (literature + arXiv survey). This document is fully
self-contained; it assumes no access to our codebase. We want ideas — algorithms,
papers, reformulations — to make the computation below dramatically faster, ideally
1–2 orders of magnitude, without sacrificing exactness.

---

## 0. One-paragraph statement of the problem

We reproduce, from first principles, the construction of a string-theory vacuum with
an exponentially small cosmological constant (McAllister–Moritz–Nally–Schachner,
arXiv:2107.09064). The final number is the vacuum energy `V0 ≈ -5.5×10⁻²⁰³` (in Planck
units; `log10|V0| ≈ -202.26`) for a specific Calabi–Yau threefold `X` with Hodge
numbers `h^{1,1}(X)=214`, `h^{2,1}(X)=4`. The expensive step is **moduli stabilization**:
solving a system of "F-flatness" equations for the 214 Kähler moduli, and then computing
worldsheet-instanton (Gopakumar–Vafa) corrections to the volume. The construction
requires choosing a **basis** of 214 divisors (out of 218) that parameterize the
superpotential. There is **no closed-form rule** for the correct basis (we proved this —
see §5); instead it is a **search**: thousands of admissible bases exist, each lands the
moduli-stabilization solution in a different "chamber" of the Kähler cone, and only some
chambers have the property that all the relevant holomorphic curves have **cheaply
computable** GV invariants. We currently **scan** admissible bases until we find such a
chamber (~rank 436 of ~4683 for this example). Each candidate costs a numerical homotopy
solve plus a large geometric coverage check. We want to go faster.

---

## 1. Mathematical / physical background (just enough to make the problem precise)

### 1.1 The geometry
- `X` is a smooth Calabi–Yau threefold realized as a hypersurface in a toric variety `V`,
  obtained from a **4-dimensional reflexive lattice polytope** `Δ°` via the
  Batyrev construction. The toric variety comes from a **fine, regular, star
  triangulation (FRST)** of `Δ°`.
- The polytope has a set of lattice points; the non-origin points that are **not interior
  to facets** become the **prime toric divisors** `D_1,…,D_n` of `X`. Here `n = h^{1,1}+4 =
  218`.
- `H^{1,1}(X) ≅ ℝ^{214}`. A choice of **divisor basis** is a choice of `h^{1,1}=214`
  of the 218 prime toric divisors that are linearly independent in `H^{1,1}` and span the
  relevant lattice. The other 4 divisors are linear combinations of the basis.
- **Triple intersection numbers** `κ_{ijk} = D_i · D_j · D_k ∈ ℤ` (a symmetric 3-tensor)
  encode the ring structure. In a 214-dim basis this is a large but sparse symmetric
  tensor.
- The **Kähler cone** `K_X ⊂ H^{1,1}` is the cone of Kähler classes; a Kähler class is
  `t = Σ t^i D_i` with all curve volumes positive. The **Mori cone** (dual) is generated
  by effective curve classes. At large `h^{1,1}` the Kähler cone has **exponentially many
  sub-chambers** (one per FRST / geometric phase); flopping a curve moves between chambers.

### 1.2 The physics we compute (per fixed flux background)
Given integer flux vectors `(K, M)` (length 4 here), an upstream "racetrack" computation
(cheap, already done) produces three scalars:
- string coupling `g_s ≈ 0.00911`,
- flux superpotential magnitude `|W0| ≈ 2.3×10⁻⁹⁰`,
- Kähler-potential factor `e^{K0} ≈ 0.2344`.

The **superpotential** is `W = W0 + Σ_i A_i exp(-c_i T_i)`, summed over a basis of 214
"rigid" divisors, where `T_i` is the complexified 4-cycle volume of divisor `i`, and
`c_i ∈ {1, 6}` (1 for a Euclidean D3-brane, 6 = dual Coxeter number of `so(8)` for a
gaugino-condensate stack on an O7-plane). The **F-flatness conditions** `∂_{T_i} W = 0`
pin each rigid divisor's 4-cycle volume to
```
    Re(T_i) = τ_i ≈ (c_i / 2π) · ln(1/|W0|).
```
Because `4-cycle volume τ_i = ½ κ_{ijk} t^j t^k` (a **quadratic** function of the 2-cycle
Kähler parameters `t`), solving for `t` is solving a **system of 214 quadratics**:
```
    ½ κ_{ijk} t^j t^k = c_i ,   i = 1..214      (the "phase-1" / naive target)
```
(then refined with α' and instanton corrections — see §1.3). This is the **core
nonlinear solve**.

The final vacuum energy is
```
    V0 = -3 · e^{K0} · (g_s⁷ / (4 V_string)²) · |W0|²,
```
where `V_string = (1/6) κ_{ijk} t^i t^j t^k − ζ(3)χ/(4(2π)³) + (GV instanton corrections)`
is the corrected Calabi–Yau volume (χ = 2(h^{1,1}−h^{2,1}) = 420). The whole point of the
paper is that with the right fluxes `|W0|` is doubly-exponentially small, so `V0` is
~10⁻²⁰³.

### 1.3 Worldsheet instanton (GV) corrections — the source of the "coverage" problem
The classical volume gets corrections from holomorphic curves (genus-0 worldsheet
instantons), weighted by **Gopakumar–Vafa invariants** `n_β` for curve classes `β`:
```
    correction ~ Σ_β n_β · Li₃(exp(-2π β·t)) / (2π)³ .
```
Only curves with **small volume** `β·t` at the solution point contribute appreciably.
"Small" = `β·t < cutoff` (cutoff = 1.0 in our runs). At the McAllister chamber there are
**344** such small curves; we must know `n_β` for each.

**Computing `n_β` is the crux.** There are three layers:
1. **Toric two-face formula** — closed-form `n_β` for curves that are toric (live in
   2-faces of the polytope). Cheap. Precomputed once into a hash table
   ("toric_gv_by_class").
2. **Minimal-degree HKTY pre-pass** — for a missed class, try its minimal generator.
3. **General GV via mirror symmetry / HKTY hypergeometric series** — exact for any class,
   but cost scales with the curve's **degree in a chosen grading**, and for high-degree
   curves this is **intractable** (we observed required degree ~700 for some chambers;
   the series would need ~700 terms in a 214-variable problem — hopeless).

**The basis/chamber determines which curves are "small," and at what degree.** A curve
`β` that is degree-1 (toric, cheap) in McAllister's chamber can be degree-700 (intractable)
in another admissible chamber. The same physical `V0` results either way (it is
basis-invariant — the stabilized point is the same geometric point), but only the "right"
chamber is **computable cheaply**. So:

> **The search is for an admissible basis whose chamber selects only small curves that
> the cheap (toric + minimal-degree) layers cover.** We call such a chamber "covered."

### 1.4 The admissibility criterion (cheap, already implemented)
The paper's §2 criterion: the basis must consist of `h^{1,1}` linearly independent **rigid**
prime toric divisors (rigidity = a combinatorial property of the divisor's lattice point;
fast to test), such that **with the basis divisor volumes set to unity, the remaining 4
divisors have strictly positive volume** — equivalently the constant covector
`τ⋆ = (1,…,1)` lies in the **dual of the effective cone** `E(X)°`. This is a cheap
`4×4` determinant / linear test per candidate. For our example there are **4683 admissible
bases**.

---

## 2. The exact computational bottleneck

For each admissible basis candidate (there are ~4683; the "covered" one is at interior-rank
~436 when candidates are ordered by how deep `τ⋆` sits in `E(X)°`), verifying whether its
chamber is "covered" currently requires:

**(A) A numerical homotopy solve of the 214-quadratic F-flat system.** We path-follow
`κ_{ijk} t^j t^k = (1−s)·(init) + s·(2 c_i)` from `s=0` to `s=1` over ~64 interpolation
steps, each step running 1–3 Newton iterations. Each Newton iteration evaluates a
`214 × 214` (basis-restricted to `214 × |basis|`) Jacobian from the sparse `κ` tensor and
solves a least-squares linear system. We run 2 initial seeds (a height-projected seed and a
scaled random seed) and select a positive-volume branch. Cost: dominant. (Plus a redundant
64-step "corrected" solve that we have since identified as **diagnostic-only** and can skip.)

**(B) A geometric "coverage" check at the solution point.** We scan **561,658 ambient
Mori-cone-cap rays** (candidate small-curve classes), compute each ray's volume `β·t`,
keep those `< cutoff`, and look each up in the toric-GV hash table; for misses we run a
minimal-degree GV pre-pass. If any selected small curve remains uncovered, the chamber is
rejected.

The shared geometric data (the 561k Mori-cap rays + toric GV table) is **basis-independent**
and computed **once** (it is itself expensive — a Double-Description-Method cone computation
over ~561k generators — but amortized).

So the per-candidate cost ≈ **(A) one ~64-step Newton homotopy on a 214-var quadratic
system + (B) a 561k-ray geometric scan with hash lookups + a minimal-degree GV pre-pass.**

We scan candidates **in interior-rank order**, stopping at the first covered chamber
(~rank 436). Wall-clock: **tens of minutes** for this single example. The production goal
is a **genetic-algorithm search over the string landscape** that calls this verification
on rare promising flux candidates (~1–2 per day), so tens of minutes per call is tolerable
but we strongly want it faster (it gates research throughput, and larger examples will be
worse).

---

## 3. Key structural facts (these constrain and enable solutions)

1. **`V0` is basis-invariant.** Every admissible covered chamber yields the same `V0`. So
   we are free to accept ANY covered chamber — we do not need McAllister's specific one.
2. **The covered/uncovered verdict is a pure function of the solution point `t` and the
   fixed (basis-independent) geometric data.** Given `t`, coverage is a cheap scan.
3. **The solution point `t` is determined by the basis** (through the target `2 c_i` and
   the chamber the homotopy lands in). Different bases → different `t` → different chamber.
4. **Coverage classification is numerically robust.** We verified (5 step counts) that the
   covered verdict is **stable** once the homotopy is converged to ≥16 of 64 steps —
   identical 344-curve set and identical `V0` at 16/32/64 steps. There were **zero** cases
   where a coarse solve falsely reported "covered." The only coarse-solve failure mode is
   the homotopy not yet reaching a positive-volume branch (then we fall back to a full
   solve — safe, just not fast).
5. **The "good vs bad" separation is stark**: a covered chamber has **0** uncovered curves;
   a bad chamber has **dozens to hundreds** (e.g. 289), at GV grading degree in the hundreds.
   There is no knife-edge.
6. **The F-flat system is a structured polynomial system**: 214 equations, each a quadratic
   form `½ κ_{ijk} t^j t^k = c_i` sharing the same sparse symmetric integer tensor `κ`. We
   seek a specific real positive-volume root in a specific cone (the Kähler cone chamber),
   not all roots.
7. **The cone/chamber structure is combinatorial**: it comes from the secondary fan of the
   polytope's triangulations; chambers are separated by flop walls (a single curve class
   going to zero volume). There may be exploitable combinatorial structure relating "which
   basis" → "which chamber" → "which curves are small."

---

## 4. What "admissible basis" and "covered chamber" mean numerically (concrete numbers)

- Example: `X` = polytope "4-214-647", `h^{1,1}=214`, `h^{2,1}=4`, χ=420.
- 218 prime toric divisors; 216 are rigid (combinatorial test); 2 are non-rigid (always
  excluded).
- Admissible bases: **4683** (rigid 214-subsets with `τ⋆ ∈ E(X)°`); enumerated by a cheap
  `4×4`-determinant test. (For other examples the count and the number of divisors to "drop"
  vary; in the worst case the drop-set is larger and the admissible count could be
  combinatorially large — a separate scaling worry.)
- McAllister's published basis drops divisors {46,130} (plus the 2 non-rigid). Its interior
  rank (by min remaining-divisor volume at `τ⋆`) is **436**.
- Shared geometric data: **561,658** ambient Mori-cone-cap rays; a toric-GV hash table.
- At the covered chamber: **344** small curves selected (all covered); `V_string ≈ 4711.4`;
  `log10|V0| = -202.26`.
- Coverage stability test (McAllister's basis, varying homotopy steps):
  | steps | result |
  |------:|--------|
  | 4, 8  | homotopy not converged (diagnostic solve fails) |
  | 16, 32, 64 | COVERED, 344 curves, `log10|V0| = -202.263` (identical) |

---

## 5. What we have tried (and the outcome of each)

### 5.1 Looking for a closed-form rule to pick the basis directly (FAILED — proven impossible)
We tried to find a cheap deterministic rule that outputs McAllister's basis (or any covered
one) without searching. We checked, across **all five** published examples
(`h^{1,1}` = 214, 113, 113, 81, 51):
- **Paper's own algorithm**: the paper (§2, "Algorithm for F-flat solutions") explicitly
  states it is a **brute-force search** over admissible bases, with the covered property
  verified **a posteriori** (instanton convergence). No closed form is given.
- **Minimum |det| of the complement** (most "unimodular" drop): FALSIFIED. Complement
  determinants across examples are {2, 6, 18, 2, −1}; McAllister ranks 307/3591/4367/11964/6
  by |det|. Not minimal.
- **"Nice" remaining volumes** (all unit / integer / balanced): FALSIFIED. One example has
  **non-integer** remaining volumes (14.33, 111.67).
- **CYTools' deterministic GLSM/HNF divisor basis**: FALSIFIED. It differs from the KKLT
  basis in all five.
- **Low/high L1-norm of the dropped lattice points**: FALSIFIED (McAllister ranks
  3359/4683 and 1096/4683 by the two norm orderings).
- **"Most interior" `τ⋆`**: FALSIFIED — matches 0 of 5 (McAllister is rank 436, not 1).

**Conclusion: there is no cheap closed-form basis rule.** The covered property is genuinely
tied to the F-flat chamber + GV coverability, which (so far) requires solving. This is the
robustly established negative result; we want either (a) a cheap predictor of coverability
we missed, or (b) a way to make the per-candidate solve+check dramatically cheaper, or (c) a
smarter search order.

### 5.2 "Any admissible basis + compute general GV" (FAILED — intractable)
Idea: since `V0` is basis-invariant, pick the most-interior admissible basis and just
compute the GV invariants of whatever curves it selects via the general HKTY/mirror-symmetry
series. FAILED: that chamber selects curves at **GV grading degree ~716**; the
hypergeometric series would need ~716 terms in 214 variables — completely intractable
(>5 min and climbing for a single curve; effectively unbounded). This is why a covered
(low-degree) chamber is required.

### 5.3 The scan (WORKS, but slow) — current baseline
Walk admissible bases in interior-rank order; for each run the full solve + coverage check;
stop at the first covered (~rank 436). Correct and reproduces `V0`. Tens of minutes.

### 5.4 Decoupling coverage from the expensive solve (PARTIAL WIN — current best)
We found that the coverage verdict needs only the **phase-1** (uncorrected) solution point,
not the full corrected pipeline, and that a **redundant** corrected "zeroth-order" homotopy
in the path was used only for a diagnostic log line (free to skip). We then built a
**two-stage scan**:
- **Stage 1 (cheap probe):** run the phase-1 homotopy at a **reduced step count** (no
  corrected solve, no intersection-tensor `χ` / B-field / target construction), then the
  cheap coverage scan. If it positively finds uncovered curves → reject the basis.
- **Stage 2 (exact):** run the full verification only on probe-survivors, in order; first
  clean one wins. A Stage-1 false-accept is caught here (fail-loud), so the result is
  exact.

Measured behavior:
- Coverage classification at reduced steps is **sound**: across 60 candidates, **0**
  false-"covered" verdicts. The probe never wrongly passes a bad chamber.
- At **24** reduced steps: ~65% of candidates cheaply rejected; ~35% had the phase-1
  homotopy not reach a positive-volume branch ("Indeterminate") → fell through to full
  verify (safe but not fast).
- At **40** reduced steps: **~94%** cheaply rejected (only ~6% Indeterminate). This is the
  current setting.

**Current status:** the two-stage scan at 40 probe-steps reproduces `log10|V0| = -202.26`
and rejects ~94% of bad chambers with a cheaper (40-step vs 64-step + corrected) solve. But
it is **still slow**, because:
1. Even a 40-step homotopy on a 214-variable quadratic system is not cheap, and we still pay
   it ~436 times (once per candidate up to the covered one).
2. The coverage scan over **561,658** rays runs per candidate, plus a minimal-degree GV
   pre-pass on the misses (which itself can be non-trivial).
3. We have not parallelized the scan across candidates at the verification level **on
   purpose** — the outer production loop (a genetic algorithm) already parallelizes across
   flux candidates / polytopes and saturates all cores, so nesting parallelism here would
   only oversubscribe. (Single-shot verification of one promising candidate, however, runs
   on essentially one core and is the slow case we care about.)

### 5.5 Things we deliberately did NOT do
- We do not "load" McAllister's basis from a data file in the search path (that would be
  cheating; the whole point is to derive everything). For exact-reproduction regression
  tests we do read it, but the GA/search must compute it.
- We do not approximate `V0` or use fallbacks; every physics quantity is computed exactly or
  the run fails loudly.

---

## 6. The specific questions we want literature/arXiv ideas on

We are looking for concrete algorithms, papers, or reformulations addressing any of:

### Q1. Cheap prediction of "covered chamber" from the basis/combinatorics.
Is there a way to predict, **without** solving the F-flat homotopy, which admissible basis
lands in a chamber whose small curves (volume `< 1` at the solution) are all toric /
low-degree? Equivalently: can we predict the **set of curves that will be small at the
F-flat point** from the basis + intersection numbers + Mori cone, cheaply? The F-flat point
is roughly where divisor volumes equal `c_i ∈ {1,6}`; is there an analytic/LP characterization
of the chamber containing the point `{½κ t t = c_i}` and of which Mori-cone rays have small
pairing there? (Relevant areas: toric geometry, secondary fans, Kähler/Mori cone chamber
structure, tropical/Newton-polytope methods, "stretched Kähler cone" literature.)

### Q2. Faster solution of the structured quadratic F-flat system.
The core solve is: find the real positive-volume root `t ∈ ℝ²¹⁴` (in a specified cone) of
`½ κ_{ijk} t^j t^k = c_i`, a system of 214 quadratics sharing one sparse symmetric integer
3-tensor. We currently use homotopy continuation (interpolating the right-hand side) with
Newton + least-squares. Questions:
- Is there a faster specialized solver for this **"diagonal-in-the-target," shared-tensor**
  quadratic system? (Numerical algebraic geometry: PHCpack, Bertini, HomotopyContinuation.jl,
  monodromy/regeneration; but we need ONE specific real root in a cone, not all complex
  roots.)
- Can the solution be **warm-started** across nearby bases? Consecutive admissible bases
  differ by swapping a few divisors; can we **continue** the solution from one basis to the
  next (a homotopy in "basis space") far cheaper than re-solving from scratch?
- Is there a good **closed-form / few-step initializer** for the positive-volume branch so we
  need far fewer Newton steps (the ~6% non-convergence at 40 steps is purely the initializer
  missing the positive branch)?
- The map `t ↦ τ_i = ½κ t t` is a quadratic "moment-like" map; is its inverse on the Kähler
  cone amenable to convex optimization (the volume `(1/6)κ ttt` is related to a convex
  function on the cone; KKLT vacua are critical points of a known functional)?

### Q3. Faster GV-coverage check.
- The per-candidate bottleneck (B) scans **561,658** Mori-cone rays. Most are toric-covered;
  only a handful of "dangerous" (non-toric, potentially high-degree) classes matter. Can we
  **precompute the set of non-toric-coverable curve classes** once, and per-candidate test
  only whether any of those is small at `t` (a tiny check) instead of scanning all 561k?
- More generally: cheap certification that "no uncovered curve is small at `t`" — e.g. an LP
  / cone-containment certificate that the dangerous classes all have `β·t ≥ cutoff`.

### Q4. Smarter search order over admissible bases.
We scan in "interior-rank" order and hit the covered basis at ~rank 436. Is there a better
ordering or a **direct construction** that reaches a covered chamber in `O(1)`/`O(log)`
candidates? E.g., is the covered chamber characterizable as the one **closest** (in some
metric) to a specific point, or reachable by a **greedy flop walk** from an easy starting
chamber toward where the `c_i`-volume point lives? (Relevant: secondary fan / flop graph
navigation, GKZ, wall-crossing.)

### Q5. Reformulating away the basis search entirely.
Is the dependence on "which 214 of 218 divisors" an artifact of the parameterization? The
physical content is: at the F-flat point, **all** rigid divisors have volume pinned to
`(c_i/2π)ln(1/W0)`; the basis is just which `h^{1,1}` we treat as independent coordinates.
Could we solve the (overdetermined-but-consistent) system over **all 216 rigid divisors**
once, directly landing on the physical point and the small-curve set, sidestepping the
per-basis chamber search? (We have not tried this; it may be the cleanest reformulation.)

### Q6. GPU / vectorization / better linear algebra.
The homotopy's inner loop is Jacobian assembly from a sparse 3-tensor + least-squares.
Are there standard accelerations (batched across candidates, GPU, randomized least-squares,
exploiting the shared sparse tensor) appropriate here?

---

## 7. Constraints any solution must respect
- **Exactness:** the final `V0` must match the published `log10|V0| = -202.26` (and the
  validated intermediates: `g_s ≈ 0.00911`, `|W0| ≈ 2.3×10⁻⁹⁰`, `V_string ≈ 4711.8`,
  `e^{K0} ≈ 0.2344`). No approximations that move these.
- **No data-file shortcuts** in the search path: everything derived from the polytope +
  fluxes.
- **Scales to larger `h^{1,1}`** (the landscape has `h^{1,1}` up to ~500+; the method must
  not blow up combinatorially — the admissible-basis count and the Mori-ray count both grow).
- Implementation is in **Rust** (we can call external solvers or port algorithms; we already
  ported the relevant CYTools/PPL pieces). Arbitrary-precision and exact integer arithmetic
  are available.

## 8. Glossary of search terms (for the literature sweep)
Calabi–Yau threefold, toric hypersurface, Batyrev mirror symmetry, reflexive polytope,
fine regular star triangulation (FRST), secondary fan, GKZ, Kähler cone, Mori cone, effective
cone, flop / flop wall / chamber, Gopakumar–Vafa invariants, HKTY (Hosono–Klemm–Theisen–Yau)
hypergeometric series, prepotential, KKLT moduli stabilization, perturbatively flat vacua
(PFV), F-flatness, racetrack, BBHL α' correction, stretched Kähler cone, CYTools, Double
Description Method, homotopy continuation, numerical algebraic geometry, polynomial system
solving, witness sets, monodromy solving, warm-start continuation, intersection numbers.
