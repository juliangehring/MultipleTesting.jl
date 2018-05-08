# MultipleTesting.jl News and Changes

## Version 0.3.0

### New Features

- Convex decreasing density estimator for π0 (#56)
- Higher criticism scores (`HigherCriticismScores`) and threshold (`HigherCriticismThreshold`) estimation (#55)
- List of publications which describe the implemented statistical methods (#54, #65)


### Changes

#### User-facing changes

- Fix p-value adjustments with unsorted inputs (#66)
- Speed up p-value combinations (#75)
- Define return types for p-value combinations (#76)
- Update documentation notebooks (#67)
- Drop support for julia 0.5, require julia 0.6 as the minimum version (#52, #64, #58)


#### Internal changes

- Simplify the method structure (#74)
- Switch to new julia syntax (#71)
- Simplify the step functions for p-value adjustments (#77)
- Update tests and references for p-value combinations (#70)
- Support the new `Distributions` interface (#52)
- Drop the `Compat` dependency (#53)
- Link comic logo in the readme file as an external resource (#69)


### Removed

- Drop `qvalues` in favour of BH-adjusted p-values (#73).


### Support

- Requires julia 0.6
- Preliminary support for julia 0.7 nightly builds


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
- Fixes Censored BUM π0 edge case when all p-values are equal to 1
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

- New π0 estimators: BUM, Censored BUM, Flat Grenander
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

- New π0 estimators: Two Step, Right Boundary
- New p-value adjustment methods: Benjamini-Hochberg Adaptive, Oracle, Sidak
- Package precompilation


### Changes / Changed

- Redefined Oracle Benjamini-Hochberg as Adaptive Benjamini-Hochberg procedure


### Support

- Requires julia 0.4
