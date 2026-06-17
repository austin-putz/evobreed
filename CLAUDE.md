# evobreed - OCS for livestock breeding

## Project Ambition

Replace common OCS software in animal breeding by creating a convenient and efficient package
to do optimal contribution selection (OCS) selection and mating plans for livestock breeders. 

## Main

This file will only be used to track important implemented functionality, if the plan 
has not been implemented, it doesn't belong in this file. 

All plans will be found in /plans/ folder within this directory. Please respect this
and don't put any future plans really within this file so we don't get mixed up later
with what is planned vs what was implemented, often there are major changes that are 
needed as we implement a new function or whatever it may be. 

## Package Versioning

Use three-part versions for releases and a fourth `.9000` component for active
development.

- `major.minor.patch`, such as `0.1.0`, is for tagged/user-facing releases.
- `major.minor.patch.9000`, such as `0.1.0.9000`, is for unreleased development
  on GitHub after the last release.
- Increment `patch` for backward-compatible bug fixes, documentation fixes that
  should be released, small maintenance updates, or dependency/check fixes.
- Increment `minor` for backward-compatible user-facing additions, new exported
  functions, meaningful workflow improvements, or larger internal changes that
  preserve existing behavior.
- Increment `major` for breaking API changes, removed or renamed exported
  functions, incompatible argument/return-value changes, or a maturity milestone
  such as the first stable `1.0.0` release.

After making a release such as `0.1.0`, immediately bump `DESCRIPTION` to the
next development version, usually `0.1.0.9000`. When preparing the next release,
drop `.9000` and choose the appropriate next release version, for example
`0.1.1` for fixes or `0.2.0` for new features.
