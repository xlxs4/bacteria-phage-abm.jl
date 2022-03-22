##
using Revise
##

##
using BenchmarkTools
using Profile
using ProfileSVG
using JET
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
        add_bacterium_single!(agent, model)
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

function p_adsorption(cell)
    nearby_phages = nearby_t(Phage, cell)
    p_host = model[cell].species === :a ? 0.9 : 0.1 # Tell the compiler that these are all ::Bacterium to help with type inference
    return p_host / (1 + exp(-length(nearby_phages)))
end

function p_lysis(phage)
    nearby_phages = nearby_t(Phage, phage)
    nearby_phages = isnothing(nearby_phages) ? 0 : length(nearby_phages)
    return 1 / (1 + model.properties.α * exp(-nearby_phages + model.properties.κ))
end

function phage_decay(phage)
    if rand() < p_phage_decay(model.decay_factor, model[phage].time_in_state)
        kill_agent!(phage, model)
    end
end

function tick_phage(phage)
    model[phage].time_in_state += 1
end

function infect(phage, cell)
    model[phage].state = :in_host
    model[phage].time_in_state = 0

    push!(model[cell].phages_inside, phage)
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
                    if rand() < model.properties.p_burst
                        println("TODO: burst!") #FIXME
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
        filter!(id -> model[id] isa t, ids)
        return ids
    end
    return single_type
end

function nearby_t(t::DataType, id)
    function keep_t_ids(pos)
        return (id -> model[id] isa t ? id : nothing).(ids_in_position(pos, model))
    end

    nearby_pos = collect(nearby_positions(model[id].pos, model, model.properties.infection_distance))
    filter!(pos -> !isempty(pos, model), nearby_pos)
    isempty(nearby_pos) && return nothing

    ids = filter(v -> !isempty(v), keep_t_ids.(nearby_pos))
    ids = reduce(vcat, ids)
    ids = convert(Vector{Int}, ids[ids.!=nothing])
    return isempty(ids) ? nothing : ids
end

function add_bacterium_single!(agent::A, model::ABM{<:Agents.DiscreteSpace,A}) where {A<:AbstractAgent}
    function random_without_bacterium(model::ABM{<:Agents.DiscreteSpace}, cutoff=0.998)
        function has_no_bacterium(pos, model)
            return !any((id -> model[id] isa Bacterium).(ids_in_position(pos, model)))
        end

        function no_bacterium_positions(model::ABM{<:Agents.DiscreteSpace})
            return Iterators.filter(i -> has_no_bacterium(i, model), positions(model))
        end

        if clamp(nagents(model) / prod(size(model.space.s)), 0.0, 1.0) < cutoff
            while true
                pos = random_position(model)
                has_no_bacterium(pos, model) && return pos
            end
        else
            no_bacterium = no_bacterium_positions(model)
            isempty(no_bacterium) && return nothing
            return rand(model.rng, collect(no_bacterium))
        end
    end

    position = random_without_bacterium(model)
    isnothing(position) && return nothing
    agent.pos = position
    add_agent_pos!(agent, model)
    return agent
end
##

##
function complex_step!(model)
    bacteria_death_inherent(by_single_type(Bacterium)(model))
    bacteria_death_lysis(by_single_type(Bacterium)(model))

    phages = by_single_type(Phage)(model)
    filter!(id -> model[id].state === :free, phages) # Tell the compiler these are all ::Phage to help with type inference
    for phage ∈ phages
        nearby_cells = nearby_t(Bacterium, phage)
        isnothing(nearby_cells) && continue
        target_cell = rand(nearby_cells)
        if rand() < p_adsorption(target_cell)
            kind = model[phage].kind
            if kind === :temperate
                if rand() < p_lysis(phage)
                    model[phage].kind = :induced_temperate
                end
            end
            if kind === :virulent || kind === :induced_temperate
                infect(phage, target_cell)
            end
        end
    end
    filter!(id -> model[id].state === :free, phages)
    tick_phage.(phages)
    phage_decay.(phages)

    model.properties.bacteria_count = length(by_single_type(Bacterium)(model))
    model.properties.phages_count = length(by_single_type(Phage)(model))
end
##

##
model = initialize()
##

##
run!(model, dummystep, complex_step!, 4;
    mdata=[:bacteria_count, :phages_count])
##