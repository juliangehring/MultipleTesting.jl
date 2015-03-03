# define the interface here..Todo mainly

abstract FdrEstimator
abstract LocalFdrEstimator <: FdrEstimator

# fit(pvalues, fdrestimator)     -> returns FdrFit object (or LocalFdrFit)

abstract FdrFit
abstract LocalFdrFit <: FdrFit

# In principle 1-1 correspondende between FdrEstimator and FdrFit types..
# E.g. GrenanderEstimator and GrenanderLocalFdrFit
# can we make use the type system to reflect this???



#FdrFit should implement

# pvalues(fdrfit) -> return p-values based on which fit was calculated
# tailfdr(fdrfit) -> return all estimated tailfdrs
# tailfdr(fdrfit, x) -> estimate tailfdr at x, where x might not have been one of the original pvals

# nulldistribution(fdrfit)
# nulldistribution(fdrfit, x)

# distribution(fdrfit)
# distribution(fdrfit,x)

# pi0(fdrfit)


# LocalFdrFit should also implement

# nulldensity
# density
# localfdr
