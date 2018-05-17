# References

## Statistical methods

### Adjustment of p-values

Arias-Castro, E., and Chen, S. (2017). Distribution-free multiple testing.
Electron. J. Statist. 11, 1983–2001.
[BarberCandes]

Barber, R.F., and Candès, E.J. (2015). Controlling the false discovery rate via
knockoffs. Ann. Statist. 43, 2055–2085.
[BarberCandes]

Benjamini, Y., and Hochberg, Y. (1995). Controlling the False Discovery Rate: A
Practical and Powerful Approach to Multiple Testing. Journal of the Royal
Statistical Society. Series B (Methodological) 57, 289–300.
[BenjaminiHochberg]

Benjamini, Y., and Liu, W. (1999). A step-down multiple hypotheses testing
procedure that controls the false discovery rate under independence. Journal of
Statistical Planning and Inference 82, 163–170.
[BenjaminiLiu]

Benjamini, Y., and Yekutieli, D. (2001). The Control of the False Discovery Rate
in Multiple Testing under Dependency. The Annals of Statistics 29, 1165–1188.
[BenjaminiYekutieli]

Bonferroni, C.E. (1936). Teoria statistica delle classi e calcolo delle
probabilita (Libreria internazionale Seeber).
[Bonferroni]

G’Sell, M.G., Wager, S., Chouldechova, A., and Tibshirani, R. (2016). Sequential
selection procedures and false discovery rate control. J. R. Stat. Soc. B 78,
423–444.
[ForwardStop]

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
significance. Biometrika 75, 800–802.
[Hochberg]

Holm, S. (1979). A Simple Sequentially Rejective Multiple Test Procedure.
Scandinavian Journal of Statistics 6, 65–70.
[Holm]

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a
modified Bonferroni test. Biometrika 75, 383–386.
[Hommel]

Šidák, Z. (1967). Rectangular Confidence Regions for the Means of Multivariate
Normal Distributions. Journal of the American Statistical Association 62,
626–633.
[Sidak]


### Estimation of π₀

Benjamini, Y., and Hochberg, Y. (2000). On the Adaptive Control of the False
Discovery Rate in Multiple Testing With Independent Statistics. Journal of
Educational and Behavioral Statistics 25, 60–83.
[LeastSlope]

Benjamini, Y., Krieger, A.M., and Yekutieli, D. (2006). Adaptive linear step-up
procedures that control the false discovery rate. Biometrika 93, 491–507.
[TwoStep]

Langaas, M., Lindqvist, B.H., and Ferkingstad, E. (2005). Estimating the
proportion of true null hypotheses, with application to DNA microarray data.
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67,
555–572.
[FlatGrenander, ConvexDecreasing]

Liang, K., and Nettleton, D. (2012). Adaptive and dynamic adaptive procedures
for false discovery rate control and estimation. Journal of the Royal
Statistical Society: Series B (Statistical Methodology) 74, 163–182.
[RightBoundary]

Markitsis, A., and Lai, Y. (2010). A censored beta mixture model for the
estimation of the proportion of non-differentially expressed genes.
Bioinformatics 26, 640–646.
[CensoredBUM]

Pounds, S., and Morris, S.W. (2003). Estimating the occurrence of false
positives and false negatives in microarray studies by approximating and
partitioning the empirical distribution of p-values. Bioinformatics 19,
1236–1242.
[BUM]

Robinson, D. (2016). Original Procedure for Choosing λ.
http://varianceexplained.org/files/pi0boot.pdf
[StoreyBootstrap]

Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control,
conservative point estimation and simultaneous conservative consistency of false
discovery rates: a unified approach. Journal of the Royal Statistical Society:
Series B (Statistical Methodology) 66, 187–205.
[Storey]


### Combination of p-values

Fisher, R.A. (1925). Statistical methods for research workers (Genesis
Publishing Pvt Ltd).
[FisherCombination]

Liptak, T. (1958). On the combination of independent tests. Magyar Tud Akad Mat
Kutato Int Kozl 3, 171–197.
[weighted StoufferCombination]

Mudholkar, G.S., and George, E.O. (1977). The Logit Statistic for Combining
Probabilities - An Overview (Rochester University NY, Dept of Statistics).
[LogitCombination]

Simes, R.J. (1986). An improved Bonferroni procedure for multiple tests of
significance. Biometrika 73, 751–754.
[SimesCombination]

Stouffer, S.A. (1949). The American soldier. Vol. 1: Adjustment during army life
(Princeton University Press).
[StoufferCombination]

Tippett, L.H.C. (1931). The Methods of Statistics. An introduction mainly for
workers in the biological sciences.
[TippettCombination]

Wilkinson, B. (1951). A statistical consideration in psychological research.
Psychological Bulletin 48, 156.
[WilkinsonCombination]


### Higher criticism

Donoho, D., and Jin, J. (2008). Higher criticism thresholding: Optimal feature
selection when useful features are rare and weak. PNAS 105, 14790–14795.
[HigherCriticismScores, HigherCriticismThreshold]

Klaus, B., and Strimmer, K. (2013). Signal identification for rare and weak
features: higher criticism or false discovery rates? Biostatistics 14, 129–143.
[HigherCriticismScores, HigherCriticismThreshold]


## Related software packages

- [PValueAdjust.jl](https://github.com/dirkschumacher/PValueAdjust.jl): P-value adjustment methods for multiple testing correction, *deprecated* [julia]

- [`p.adjust`](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/p.adjust.html) in the [stats](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/00Index.html) package [R]

- [mutoss](https://cran.r-project.org/web/packages/mutoss/index.html): Multiple hypothesis testing procedures [R]

- [qvalue](https://bioconductor.org/packages/release/bioc/html/qvalue.html): Q-value estimation for false discovery rate control [R]

- [fdrtool](https://cran.r-project.org/web/packages/fdrtool/index.html): Estimation of (local) false discovery rates and higher criticism [R]

- [pi0](https://cran.r-project.org/web/packages/pi0/index.html): Estimating the proportion of true null hypotheses for FDR, *deprecated* [R]
