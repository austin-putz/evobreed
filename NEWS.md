# evobreed 0.1.0.9000

## Development version

- Added the initial R package skeleton for `evobreed`.
- Added exported pedigree utilities:
  - `build_ped()` for multi-generation pedigree simulation.
  - `build_a_matrix()` for additive relationship matrix construction, using `pedigreemm` when available and a Henderson tabular-method fallback otherwise.
- Moved optimization example workflows to `scripts/`.
- Added package documentation, tests, README installation instructions, and GPL-3 licensing.

## Versioning

evobreed uses three-part versions for releases and a fourth `.9000` component
for unreleased development versions.

- Release versions use `major.minor.patch`, for example `0.1.0`.
- Development versions use `major.minor.patch.9000`, for example
  `0.1.0.9000`.
- Patch releases increment the third component, for example `0.1.1`.
- Feature releases increment the second component, for example `0.2.0`.
- Breaking or maturity-signaling releases increment the first component, for
  example `1.0.0`.

The `.9000` convention is common in R package development because it makes
GitHub/development installs sort after the last release but before the next
release. For example, after releasing `0.1.0`, development work continues as
`0.1.0.9000`; the next public release would usually become `0.1.1` for fixes
or `0.2.0` for new user-facing features.
