# Proven Difficulty Lever Implementations

Pick one or more per subproblem in Phase 1 Step 2. These are patterns, not a checklist, and not mutually exclusive — accepted tasks regularly stack two or three of these in the same section (a boundary trap, a complex-phase trap, and an aggregate-vs-count trap have all appeared together in one data-processing subproblem). Adapt the specifics to the paper's actual physics rather than forcing a lever that doesn't fit naturally.

## Lever A: Multi-part return with mixed types (SP1)

```python
# Prompt says: return (J, weight_folded_low, weight_folded_high, n_folded_low, n_folded_high)
# weights are float (spectral sums), counts are int (mode counts)
# Discriminates: model returns a shorter tuple, swaps weight/count positions, or gets the types wrong
```

## Lever B: n-vs-n-2 / missing-factor denominator traps (SP2)

```python
# Variant 1 — a separate weighted_rmse output:
# Add to the return dict: 'weighted_rmse': float(np.sqrt(2.0 * result.cost / n_points))
# Trap: the natural implementation computes sqrt(chi2_red) = sqrt(2*cost/(n-2)) instead of sqrt(2*cost/n)
# With n=40 this is wrong by 2.6%, which rtol=1e-2 catches cleanly
# Discriminative test:
#   expected = np.sqrt(2.0 * out['cost'] / n)
#   wrong = np.sqrt(out['chi2_red'])
#   assert np.isclose(out['weighted_rmse'], expected, rtol=1e-3)
#   assert not np.isclose(out['weighted_rmse'], wrong, rtol=1e-3)

# Variant 2 — chi2_red's own missing factor of 2:
# chi2_red = 2.0 * result.cost / dof            (dof = n_points - 2, since least_squares'
# .cost is already one-half the sum of squared residuals)
# Trap: the natural implementation forgets the factor of 2 and reports cost/dof instead
# Discriminative test:
#   expected = 2.0 * out['cost'] / (n - 2)
#   wrong = out['cost'] / (n - 2)
#   assert np.isclose(out['chi2_red'], expected, rtol=1e-9), 'chi2_red must equal 2*cost/(n-2)'
#   assert not np.isclose(out['chi2_red'], wrong, rtol=1e-3), 'factor-of-2 must be present'
```

## Lever C: Zero-cost correlation must stay finite (SP2)

```python
# rho_ae_ah must come from the UNSCALED (J^T J)^-1, not the chi2-scaled covariance
# At zero cost: chi2=0, the scaled covariance is all zeros, so the naive rho = 0/0 = NaN
# Correct approach: use the unscaled shape matrix -> rho ~ -0.335 (finite, negative)
# Discriminative test:
#   assert np.isfinite(out['rho_ae_ah'])
#   assert out['rho_ae_ah'] < -0.1  # must be negative, not 0.0 or NaN
```

## Lever D: None handling — catch both exception types (SP1/SP2)

```python
# The test must accept EITHER ValueError OR TypeError:
#   try:
#       my_function(..., param=None)
#   except (ValueError, TypeError):
#       failed = True
# Why: the natural implementation does np.isfinite(None) -> TypeError, not the spec's ValueError
# If the test only catches ValueError, the model fails this 100% of the time -> too hard
# If it catches both, the model only fails on the other discriminators -> stays in band
```

## Lever E: Self-contained MP (MP)

```python
# The MP SOLUTION must not call SP1/SP2 as external imports — it defines working copies
# of them as its own module-level helpers and calls those (per Phase 2's Solution rules).
# This forces the model to embed correct copies of SP1 + SP2 from memory.

# The MP TESTS follow a matching rule on the *test* side: don't call SP1/SP2 inside the
# test either. In practice this means one of two things:
#   (a) hand-derive the expected value inline from the raw inputs, e.g.:
#         g_l = cp * (DC * re_l - DV * rh_l)
#         g_r = cp * (DC * re_r - DV * rh_r)
#         wgt = np.abs(g_l - g_r) ** 2
#         idx = np.clip(np.searchsorted(omega_grid, omega_k, 'right') - 1, 0, n_bins - 1)
#         J_num = np.zeros(n_bins); np.add.at(J_num, idx, wgt); J_num /= np.diff(omega_grid)
#       — rare in practice, since a full MP output is usually too composite for this; or
#   (b) hardcode a high-precision expected value obtained by running the correct reference
#       solution once, and assert against it with a tight tolerance, e.g.:
#         assert np.isclose(out['enhancement_factor'], 1.4467501373897311, rtol=1e-7, atol=1e-12)
#       (b) is the more common pattern in accepted tasks — a golden-value regression test,
#       not a from-scratch recomputation. Use a fixed RNG seed for the input data so the
#       golden value is reproducible.
```

## Lever F: Complex-valued quantities where the phase carries the trap (SP1 / data processing)

```python
# When a coupling or amplitude is complex (dtype complex128), a naive implementation
# often works only with the real part, or takes the magnitude of each term separately
# before combining them, instead of the squared magnitude of the combined complex
# difference.
# Discriminative test: give the inputs nonzero imaginary parts and check the squared
# magnitude of the complex difference against the real-part-only shortcut:
#   g = cp * (Dc * rho_lambda - Dv * rho_ref)          # complex
#   expected = abs(g) ** 2
#   wrong = (cp * (Dc * rho_lambda.real - Dv * rho_ref.real)) ** 2
#   assert np.isclose(J[0], expected, rtol=1e-10)
#   assert not np.isclose(J[0], wrong, rtol=1e-3), 'imaginary parts must contribute'
```

## Lever G: Boundary-inclusion convention (binning/histogram tasks)

```python
# A half-open binning convention ([lo, hi) per bin) is standard, but the LAST bin is
# often closed at both ends so a value exactly on the outer edge is still counted,
# rather than falling out of range. A naive implementation applies the half-open rule
# uniformly and drops (or wrongly folds) a mode landing exactly on the final edge.
# Discriminative test: place a mode exactly at omega_grid[-1] and check it lands in
# the last bin as an in-range value, not a folded one:
#   assert np.isclose(J[-1], expected_value_including_that_mode, rtol=1e-10), \
#       'last bin must be closed at both ends'
#   assert wfh == 0.0 and nh == 0, 'mode exactly on the edge is in range, not folded'
```

## Lever H: Aggregate weight vs. bare count confusion (diagnostic outputs)

```python
# When a function reports both a count (e.g. n_folded_low, an int) and a physically
# weighted sum over the same subset (e.g. weight_folded_low, a float sum of |g|^2), a
# naive implementation sometimes conflates the two -- e.g. returning the count cast to
# float instead of the actual weighted sum.
# Discriminative test:
#   expected_weight = sum of |g_k|^2 over the folded modes  # NOT just their count
#   assert np.isclose(wfl, expected_weight, rtol=1e-10)
#   assert not np.isclose(wfl, float(n_folded_low), rtol=1e-3), \
#       'weight_folded_low is a sum of |g|^2, not a mode count'
```
