
# Comparing estimators for π0

## Setting the scene


```julia
using MultipleTesting
```


```julia
using Distributions
using DataFrames
using Gadfly
using Compose
```

## Simulation

We will simulate p-values under the null and alternative hypothesis by drawing from random variable distributions. For the null hypothesis, p-values are uniformly distributed in the interval [0,1]; p-values under the alternative hypothesis are modelled by a Beta distribution with a momotonically decreasing density.


```julia
π0 = 0.7
bum = BetaUniformMixtureModel(π0, 0.5, 3.0)
```




    MixtureModel{Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}}(K = 2)
    components[1] (prior = 0.3000): Distributions.Beta(α=0.5, β=3.0)
    components[2] (prior = 0.7000): Distributions.Uniform(a=0.0, b=1.0)





```julia
h1, h0 = components(bum);
```


```julia
x = linspace(0.01, 1.0, 200)
```




    linspace(0.01,1.0,200)




```julia
p = plot(
    layer(x = x, y = pdf(bum, x), Geom.line),
    layer(x = x, y = pdf(h0, x), Geom.line, Theme(default_color = colorant"darkgreen")),
    layer(x = x, y = pdf(h1, x), Geom.line, Theme(default_color = colorant"darkred")),
    Guide.xlabel("p-value"), Guide.ylabel("density")
)
draw(PNG(20cm, 10cm), p)
```


![png](pi0-estimation_files/pi0-estimation_9_0.png)


## Assessment of π0 estimators

We assess the performance of three estimators for π0, the fraction of tests under the null hypothesis. Here, we compare different estimators. Using the mixture model from above, we simulate p-values with a true π0 of 0.8, i.e. 20% of the p-values come from the alternative model.


```julia
estimators = ["Storey",
              "StoreyBootstrap",
              "LeastSlope",
              "TwoStep",
              "RightBoundary",
              "CensoredBUM",
              "BUM",
              "FlatGrenander"];
```


```julia
m = 500
```




    500




```julia
pi0hat = zeros(m, length(estimators))
for i in 1:m
    pvals = rand(bum, m)
    for (j, e) in enumerate(estimators)
        pi0estimator = call(eval(parse(e)))
        pi0hat[i,j] = estimate_pi0(pvals, pi0estimator)
    end
end
```


```julia
df = convert(DataFrame, pi0hat)
names!(df, map(symbol, estimators))
df = stack(df, collect(1:length(estimators)))
names!(df, [:estimator, :pi0hat]);
```


```julia
p = plot(
    layer(yintercept = [π0], Geom.hline(color = colorant"black")),
    layer(df, x = "estimator", y = "pi0hat", color = "estimator", Geom.boxplot)
)
draw(PNG(20cm, 10cm), p)
```

Visualizing the estimated π0 for the three estimators shows us a typical bias-variance tradeoff: The bootstrap has a small bias with regard to the true value (black horizontal line), but exhibits a large degree of variablity. In contrast, the two step estimator has the lowest spread, while clearly overestimating π0.
