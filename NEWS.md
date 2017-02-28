# MultipleTesting.jl News and Changes

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


### Changes

- Optimised Benjamini-Liu π0 estimator and test cases
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
- New p-value adjustment methods: Forward Step, Benjamini-Liu
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
