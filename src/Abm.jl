module Abm

# TODO: make this cleaner, see https://github.com/sbromberger/LightGraphs.jl/blob/master/src/LightGraphs.jl
using Agents:
    @agent, ABM, AbstractAgent, DiscreteSpace, GridSpace, add_agent!,
    add_agent_pos!, allids, dummystep, genocide!, ids_in_position, kill_agent!,
    move_agent!, nagents, nearby_ids, nearby_positions, positions, random_position, run!

using InteractiveDynamics: abmplot
using GLMakie
using Random: MersenneTwister, rand

include("impl.jl")
include("parameters.jl")
include("step.jl")
include("utils.jl")

function initialize(; M=33, seed=125)
    space = GridSpace((M, M))
    properties = Parameters()
    rng = MersenneTwister(seed)

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
    params = Dict(
        :a => 0.0:0.1:1.0,
        :b => 0.0:0.1:1.0,
        :m => 0.0:1.0:40.0,
        :α => 0.0:1.0:100.0,
        :κ => 0.01:0.01:0.5,
        :moi_proxy_radius => 1:3,
        :infection_distance => 1:3,
        :latent_period => 1:10,
        :burst_size => 3:8,
        :growth_rate => 0.2:0.1:1.0,
        :decay_factor => 0.01:0.01:0.2,
        :p_burst => 0.4:0.1:0.9
    )

    ac(a::A) where {A<:AbstractAgent} = a.type === :bacterium ? :cyan : :magenta
    as = 10
    am(a::A) where {A<:AbstractAgent} = a.type === :bacterium ? :circle : :diamond
    f(model::ABM{<:DiscreteSpace,A}) where {A<:AbstractAgent} = map(v -> length(v), model.space.s)

    heatarray = f
    heatkwargs = (colorrange=(0, 7), colormap=:thermal)
    plotkwargs = (;
        ac, as, am,
        scatterkwargs=(strokewidth=1.0,),
        heatarray, heatkwargs
    )

    fig, ax, abmobs = abmplot(model; (model_step!)=complex_step!, params, plotkwargs...)
    return fig
end

function run(model, n)
    run!(model, dummystep, complex_step!, n;
        mdata=[:bacteria_count, :phages_count])
end

function main()
    model = initialize()
    run(model, 399)
    plot(model)
end

end