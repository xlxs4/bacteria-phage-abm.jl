module Abm

# TODO: make this cleaner, see https://github.com/sbromberger/LightGraphs.jl/blob/master/src/LightGraphs.jl
using Agents, InteractiveDynamics, Random, CairoMakie

include("impl.jl")
include("parameters.jl")
include("step.jl")
include("utils.jl")

function initialize(; M=33, seed=125)
    space = GridSpace((M, M))
    properties = Parameters()
    rng = Random.MersenneTwister(seed)

    model = ABM(
        Organism, space;
        properties, rng
    )

    for n ∈ 1:properties.bacteria_count
        agent = Organism(n, (1, 1), rand(model.rng) < 0.5 ? :a : :b, Vector{Int}(), :nothing, :nothing, -1, :bacterium)
        add_bacterium_single!(agent, model)
    end
    for _ ∈ 1:properties.phages_count
        roll = rand(model.rng)
        if roll < 0.45
            kind = :temperate
        elseif roll < 0.85
            kind = :virulent
        else
            kind = :deficient
        end
        add_agent!(Organism, model, roll < 0.5 ? :a : :b, Vector{Int}(), kind, :free, 0, :phage)
    end

    return model
end

function plot(model::ABM)
    ac(a::A) where {A<:AbstractAgent} = a.type === :bacterium ? :cyan : :magenta
    as = 10
    am(a::A) where {A<:AbstractAgent} = a.type === :bacterium ? :circle : :diamond
    f(model::ABM{<:Agents.DiscreteSpace,A}) where {A<:AbstractAgent} = map(v -> length(v), model.space.s)

    heatarray = f
    heatkwargs = (colorrange=(0, 7), colormap=:thermal)
    plotkwargs = (;
        ac, as, am,
        scatterkwargs=(strokewidth=1.0,),
        heatarray, heatkwargs
    )

    fig, ax, abmobs = abmplot(model; model_step! = complex_step!, plotkwargs...)
    return fig
end

function main()
    model = initialize()

    run!(model, dummystep, complex_step!, 399;
        mdata=[:bacteria_count, :phages_count])

    plot(model)
end

end