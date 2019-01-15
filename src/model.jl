"""
Beta Uniform Mixture (BUM) Model


# Arguments

- π0 : Contributing fraction of the uniform distribution to the full model
- α, β : Parameters of the Beta distribution, Float64, default: 0.5, 3.0


# Return values

`MixtureModel`, as defined in the `Distributions` package, composed of

- a uniform distribution in the interval [0, 1], with weight/prior π₀
- a Beta distribution with parameters α and β, with weight/prior 1-π₀


# Examples

```jldoctest
julia> bum = BetaUniformMixtureModel(0.2, 0.5, 1.0);

julia> using Distributions

julia> pdf.(bum, 0.2:0.2:1.0)
5-element Array{Float64,1}:
 1.094427190999916
 0.832455532033676
 0.7163977794943224
 0.647213595499958
 0.6000000000000001

```


# References

Pounds, S., and Morris, S.W. (2003). Estimating the occurrence of false
positives and false negatives in microarray studies by approximating and
partitioning the empirical distribution of p-values. Bioinformatics 19,
1236–1242.

"""
function BetaUniformMixtureModel(π0::AbstractFloat, α::AbstractFloat = 0.5, β::AbstractFloat = 1.0)
    isin(π0, 0, 1) || throw(DomainError("π0 must be in [0, 1]"))
    bum = MixtureModel([Beta(α, β), Uniform()], [1-π0, π0])
    return bum
end
