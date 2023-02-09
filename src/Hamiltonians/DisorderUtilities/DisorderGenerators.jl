
using Random
using Distributions

export GaussianDisorder, BoxDisorder, QuenchedDisorderGenerator, generate_disorder

abstract type DisorderDistribution end
struct GaussianDisorder <: DisorderDistribution end
struct BoxDisorder <: DisorderDistribution end

distribution(dis::DisorderDistribution, args...) = throw(MethodError(distribution, (dis, args...)))
distribution(dis::GaussianDisorder, args...) = Distributions.Normal(args...)
distribution(dis::BoxDisorder, center, width) = Distributions.Uniform(center - width/2, center + width/2)

struct QuenchedDisorderGenerator{T <: AbstractFloat, D <: DisorderDistribution, G <: Random.AbstractRNG}
    center::T
    width::T
    rng::G
    disorder_distribution::D
end

function QuenchedDisorderGenerator(rng::Random.AbstractRNG, center, width, dis::DisorderDistribution)
    tup = promote(center, width)
    return QuenchedDisorderGenerator{eltype(tup), typeof(dis), typeof(rng)}(tup..., rng, dis)
end
QuenchedDisorderGenerator(center, width, dis::DisorderDistribution) = QuenchedDisorderGenerator(Random.GLOBAL_RNG, center, width, dis)

generate_disorder( qdg::QuenchedDisorderGenerator, num_or_size ) = rand(qdg.rng, distribution(qdg.disorder_distribution, qdg.center, qdg.width), num_or_size)