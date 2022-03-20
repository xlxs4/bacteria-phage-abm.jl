##
using Revise
##

##
using Agents
using Random
##

##
@agent Bacterium GridAgent{2} begin
    species::Symbol
    phages_inside::Vector{Int}
end

@agent Phage GridAgent{2} begin
    species::Symbol
    kind::Symbol
    state::Symbol
    time_in_state::Int
end

Base.@kwdef mutable struct Parameters
    bacteria_count::Int = 265
    phages_count::Int = 148
    environment::Symbol = :semi_solid
    diffuse::Bool = true
    a::Float64 = 0.0
    b::Float64 = 0.08
    m::Float64 = 27.0
    α::Float64 = 8
    κ::Float64 = 0.05
    moi_proxy_radius::Int = 1
    infection_distance::Int = 1
    latent_period::Int = 5
    burst_size::Int = 4
    carrying_capacity::Int = 3 * bacteria_count
    growth_rate::Float64 = 0.6
    decay_factor::Float64 = 0.2
    p_burst = 0.9
end
##

##
function initialize(; N=200, M=20, seed=125)
    space = GridSpace((M, M))
    properties = Parameters()
    rng = Random.MersenneTwister(seed)

    model = ABM(
        Union{Bacterium,Phage}, space;
        properties, rng
    )

    for n ∈ 1:N/2
        roll = rand()
        agent = Bacterium(n, (1, 1), roll < 0.5 ? :a : :b, Vector{Int}())
        add_agent_single!(agent, model)
    end
    for n ∈ N/2+1:N
        roll = rand()
        if roll < 0.45
            kind = :temperate
        elseif roll < 0.85
            kind = :virulent
        else
            kind = :deficient
        end
        add_agent!(Phage, model, roll < 0.5 ? :a : :b, kind, :free, 0.0)
    end

    return model
end
##

##
p_death(a, b, m) = a + ((1 - a) / (1 + exp(-b * (-m))))

p_phage_decay(decay_factor, time_in_state) = 1 - (1 * exp(-decay_factor * time_in_state))

function phage_decay(phage)
    if rand() < p_phage_decay(model.decay_factor, model[phage].time_in_state)
        kill_agent!(phage, model)
    end
end

function tick_phage(phage)
    model[phage].time_in_state += 1
end

function bacteria_death_inherent(bacteria)
    for id ∈ bacteria
        if rand() < p_death(model.a, model.b, model.m)
            genocide!(model, model[id].phages_inside)
            kill_agent!(id, model)
        end
    end
end

function bacteria_death_lysis(bacteria)
    for id ∈ bacteria
        phages_inside = model[id].phages_inside
        if !isempty(phages_inside)
            tick_phage.(phages_inside)
            for phage ∈ phages_inside
                agent = model[phage]
                if (agent.time_in_state ≥ model.latent_period) && (agent.kind === :virulent || agent.kind === :induced_temperate)
                    if rand() < p_burst
                        println("TODO: burst!")
                    end
                end
            end

            isempty(model[id].phages_inside) && kill_agent!(id, model)
        end
    end
end
##

##
function by_single_type(t::DataType)
    function single_type(model::ABM)
        ids = collect(allids(model))
        filter!(id -> typeof(model[id]) == t, ids)
        return ids
    end
    return single_type
end
##

##
function complex_step!(model)
    bacteria_death_inherent(by_single_type(Bacterium)(model))
    bacteria_death_lysis(by_single_type(Bacterium)(model))

    ids = by_single_type(Phage)(model)
    filter!(id -> model[id].state === :free, ids)
    phage_decay.(ids)

    model.properties.bacteria_count = length(by_single_type(Bacterium)(model))
    model.properties.phages_count = length(by_single_type(Phage)(model))
end
##

##
model = initialize()
##

##
run!(model, dummystep, complex_step!, 9; mdata=[:bacteria_count])
##