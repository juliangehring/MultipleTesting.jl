# MultipleTesting.jl News and Changes

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
