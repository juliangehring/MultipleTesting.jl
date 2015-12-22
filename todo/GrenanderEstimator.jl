type GrenanderEstimator <: LocalFdrEstimator
  pi0estimator:: Pi0Estimator
end



# saves everything in sorted way for now, change that later...
type GrenanderLocalFdrFit <: LocalFdrFit
  pvalues   :: Array{Float64,1}
  f         :: Array{Float64,1} #density
  F         :: Array{Float64,1} #distribution
  pi0       :: Float64
end

pvalues(fdrfit::GrenanderLocalFdrFit) = fdrfit.pvalues
pi0(fdrfit::GrenanderLocalFdrFit)  = fdrfit.pi0
distribution(fdrfit::GrenanderLocalFdrFit)    = fdrfit.F
density(fdrfit::GrenanderLocalFdrFit)    =    fdrfit.f
localfdr(fdrfit::GrenanderLocalFdrFit)   =    fdrfit.pi0./fdrfit.f #TODO make sure between [0,1]
tailfdr(fdrfit::GrenanderLocalFdrFit)    =    fdrfit.pi0.*fdrfit.pvalues./fdrfit.F #TODO make sure between [0,1]


function fit(pv:: Array{Float64,1}, fdrestimator::GrenanderEstimator)
  m = length(pv)
  pi0 = estimatepi0(pv, fdrestimator.pi0estimator)

  pv_ecdf_fun = ecdf(pv)          #going about it in quite inefficient way but this should not be the bottleneck anyway
  pv_sorted = sort(pv)
  ecdf_value = Array(Float64,m)

  # apply modification to make grenander estimate consistent with two-groups model
  for i=1:m
    ecdf_value[i] = pv_ecdf_fun(pv_sorted[i])
    if ecdf_value[i] < pi0*pv_sorted[i]
      ecdf_value[i] = pi0*pv_sorted[i]
    elseif ecdf_value[i] > 1-pi0*(1-pv_sorted[i])
      ecdf_value[i] = 1-pi0*(1-pv_sorted[i])
    end
  end

  #denote by (x,y) the coords of the modified ecdf
  Δx = diff(pv_sorted)
  Δy = float64(diff(ecdf_value))

  f = Δy./Δx #rough density estimate


  f = -isotonicregression(-f, Δx) # make previous density estimate obey monotonicity
  F  = ecdf_value[1] + vcat(0,cumsum(f .* Δx))
  f = push!(f, f[m-1])


  GrenanderLocalFdrFit(pv_sorted,f,F,pi0)
end
