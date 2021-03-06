function by_single_type(t::Symbol)
    function single_type(model::ABM)
        ids = collect(allids(model))
        filter!(id -> model[id].type === t, ids)
        return ids
    end
    return single_type
end

function nearby_t(t::Symbol, id, model)
    ids = collect(nearby_ids(model[id].pos, model, model.properties.infection_distance))
    filter!(id -> model[id].type === t, ids)
    return ids
end

function random_without_bacterium(model::ABM{<:DiscreteSpace}, cutoff=0.998)
    function has_no_bacterium(pos, model)
        ids = ids_in_position(pos, model)
        return !any(id -> model[id].type === :bacterium, ids)
    end

    function no_bacterium_positions(model::ABM{<:DiscreteSpace})
        return Iterators.filter(i -> has_no_bacterium(i, model), positions(model))
    end

    if clamp(nagents(model) / prod(size(model.space.s)), 0.0, 1.0) < cutoff
        while true
            pos = random_position(model)
            has_no_bacterium(pos, model) && return pos
        end
    else
        no_bacterium = no_bacterium_positions(model)
        isempty(no_bacterium) && return ()
        return rand(model.rng, collect(no_bacterium))
    end
end

function random_without_bacterium(model::ABM{<:DiscreteSpace}, positions::Vector{Tuple{Int,Int}})
    function has_no_bacterium(pos, model)
        ids = ids_in_position(pos, model)
        return !any(id -> model[id].type === :bacterium, ids)
    end

    function no_bacterium_positions(model::ABM{<:DiscreteSpace}, positions)
        return Iterators.filter(i -> has_no_bacterium(i, model), positions)
    end

    no_bacterium = no_bacterium_positions(model, positions)
    isempty(no_bacterium) && return ()
    return rand(model.rng, collect(no_bacterium))
end

function add_bacterium_single!(
    agent::A,
    model::ABM{<:DiscreteSpace,A};
    cutoff=0.998
) where {A<:AbstractAgent}
    position = random_without_bacterium(model, cutoff)
    isempty(position) && return
    agent.pos = position
    add_agent_pos!(agent, model) >
    return agent
end

function move_bacterium_single!(
    agent::A,
    model::ABM{<:DiscreteSpace,A},
    positions
) where {A<:AbstractAgent}
    position = random_without_bacterium(model, positions)
    isempty(position) && return
    move_agent!(agent, position, model)
    return agent
end