# Upstream integration — November 15, 2026

Merged source: scottprahl/iadpython main at `33885de93c5d991c518213951ec4890011dc90e0`, into the fork based on published master `300e3029fb663fabb3d15b3e130ee9358e36e1cc`.

The preparation commit normalizes Python line endings, declares pandas, and restores the legacy HG redistribution entry point. The merge preserves the fork's Layer API, mismatched-index transport, tabulated phase functions, moment inputs and custom delta-M calculations while importing upstream inverse, sphere, MC, CLI and packaging work.

## Narrow integration changes

- Added upstream's forward-call counter without replacing the fork's scalar, spectral or layered calculation paths.
- New upstream measurement helpers use the Sample ballistic solver, preserving the fork's external-angle convention and full-stack transport. Flat ballistic calculations now pass slide optical depths explicitly.
- Default boundary index follows updates to the transport index. An explicitly supplied boundary index remains independent. This fixes upstream loaders that change `sample.n` after construction.
- `Experiment.point_at` slices phase tables, moments and layer wavelength data and clears phase caches.
- Preserved `verbose` compatibility alongside upstream hot-start and simplex arguments.
- Fixed exactly matched Fresnel interfaces at grazing incidence. Updated two old specular tests to pass external incidence angles; their independent physical reference values are unchanged.
- Added nbformat as a runtime dependency for upstream's iadplus command; preserved all new entry points and packaged scripts.
- Replaced upstream's Unix-only date formatting with equivalent portable formatting for Windows CLI output.

The fork's semi-infinite coefficient normalization remains unchanged intentionally. It differs from upstream's scattering-normalized convention. No claim is made that MC corrections or inversion support arbitrary custom-phase, mismatched multilayer specimens.

## Preservation evidence

The focused suite initially passed 149 tests and two subtests. Subsequent integration tests cover default-boundary synchronization and frozen pre-merge reference results; all 11 integration tests pass.

The frozen reference in `tests/data/fork_forward_reference.json` was generated from the published pre-merge fork, not the modified implementation. It covers three-layer stacks mixing HG, TABULATED and MOMENTS, two wavelengths, unequal indices, and both air and asymmetric external media. All four R/T outputs and ballistic R/T matched with maximum absolute difference **0.0** in the direct comparison. The committed test allows small floating-point differences across environments.

The wheel built successfully and was installed in an isolated uv environment outside the checkout. Imports of iadpython, pandas and nbformat succeeded, all three console entry points were present, and an unequal-index layered solve ran successfully from the installed wheel.

The previously reported lossless-stack energy deficit remains a known numerical limitation; this merge does not attempt to redesign that solver.

## Final validation

The broad non-slow run completed in 14 minutes: 402 passed, 28 skipped, 43 subtests passed, with six failures and five setup errors. Those remaining cases were investigated and resolved as follows:

- Two CLI tests: use platform-aware quoted path expectations, and explicitly disable MC in the pure inversion smoke test. Both pass on rerun.
- Five C-parity setup errors: the upstream library loader exits before its intended skip check when libiad is unavailable. Test setup now reports a skip in that situation; it does not alter scientific code.
- Four data2 tests: upstream does not include the optional tests/data2 directory in Git. These classes now skip when that directory is absent; they still execute normally when it is supplied.

The final targeted rerun returned **13 passed, 9 skipped**. This includes all affected tests and all 11 integration tests, including one frozen-reference test added after the broad run was collected.

Across the broad run and targeted reruns, final coverage is **405 distinct tests passed, 37 skipped, and 43 subtests passed**. One slow test was deselected, and notebooks were excluded by the repository's existing pytest configuration. There are no unresolved failures among the executed checks. External C/MC computations and unavailable optional fixtures are not claimed as validated.

The broad command was `python -u -m pytest tests -v -m "not slow" --timeout=120 --tb=short`, in an isolated Python 3.12 uv environment with numpy, scipy, pandas, matplotlib, pytest, pytest-timeout and nbformat. A targeted rerun covered test_merge_compatibility.py, both iadc parity files, test_rxt2.py, and the two corrected CLI tests.

The upstream commit was rechecked before finalization and remained 33885de. The merge retains both parents so Git records all upstream history as incorporated; unrelated local untracked files are excluded from the commits.

