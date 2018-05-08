# MultipleTesting.jl News and Changes

## Unreleased

### Changes

- Drop support for julia 0.5, require julia 0.6 as minimum version
- Drop `Compat` dependency
- Support new `Distributions` interface


### Support

- Requires julia 0.6


## Version 0.2.3

### Changes

- Switch to new `Weights` type in 'StatsBase' (>= v0.15.0)


## Version 0.2.2

### Changes

- Makes the `PValues` type immutable and persistent
- Finalises support of julia 0.6 with testing


### Support

- Requires julia 0.5 or 0.6


## Version 0.2.1

### Changes

- Adds compatibility with julia 0.6 development versions, including automated tests
- Fixes Censored BUM π₀ edge case when all p-values are equal to 1
- Updates test cases and reference values for `qValues`
- Introduces `Compat` as a new dependency
- Separates test coverage reports for different julia versions


### Support

- Requires julia 0.5 or 0.6


## Version 0.2.0

### New Features

- Combination of p-values through the `combine` interface:
  * Fisher
  * Stouffer
  * Logit
  * Wilkinson
  * Simes
  * Tippett
  * Minimum of adjusted p-values
- Analysis notebook for p-value combinations
- `PValues` type: Representation of vectors with p-values
- Adjustments of p-values support specifying the total number of tests
- Barber-Candès p-value adjustment


### Changes

- Optimised Benjamini-Liu p-value adjustment and test cases
- Estimator types are now immutables
- Harmonised type names
- Restructed readme
- Simplified method signatures
- Various performance improvements


### Removed

- Export of lower-level functions for p-value adjustment


### Support

- Requires julia 0.5


## Version 0.1.0

### New Features

- New π₀ estimators: BUM, Censored BUM, Flat Grenander
- New p-value adjustment methods: Forward Stop, Benjamini-Liu
- BUM model for simulating p-value distributions
- Fast isotonic regression
- Reference list of related software packages in documentation


### Changes

- Grenander ECDF estimator
- Immutable method types
- Default α value for Two Step constructor
- Tests use new julia 0.5 base testing framework
- Documentation moved to 'docs' directory


### Removed

- Test coverage reporting through coveralls


### Support

- Requires julia 0.5


## Version 0.0.2


### New features

- New π₀ estimators: Two Step, Right Boundary
- New p-value adjustment methods: Benjamini-Hochberg Adaptive, Oracle, Sidak
- Package precompilation


### Changes / Changed

- Redefined Oracle Benjamini-Hochberg as Adaptive Benjamini-Hochberg procedure


### Support

- Requires julia 0.4
