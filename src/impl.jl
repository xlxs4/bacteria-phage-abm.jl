p_death(a, b, m) = a + ((1 - a) / (1 + exp(-b * (-m))))

p_phage_decay(decay_factor, time_in_state) = 1 - (1 * exp(-decay_factor * time_in_state))

function p_lysis(phage, model)
    nearby_phages = nearby_t(:phage, phage, model)
    nearby_phages = length(nearby_phages)
    return 1 / (1 + model.properties.α * exp(-nearby_phages + model.properties.κ))
end

function p_adsorption(cell, model)
    nearby_phages = nearby_t(:phage, cell, model)
    p_host = model[cell].species === :a ? 0.9 : 0.1
    return p_host / (1 + exp(-length(nearby_phages)))
end

function p_grow(model)
    cells = by_single_type(:bacterium)(model)
    filter!(id -> isempty(model[id].phages_inside), cells)
    cell_count = length(cells)

    properties = model.properties
    ΔN = properties.growth_rate * cell_count * (1 - cell_count / properties.carrying_capacity)
    return ΔN / cell_count
end

function phage_decay(phage, model)
    if rand(model.rng) < p_phage_decay(model.decay_factor, model[phage].time_in_state)
        kill_agent!(phage, model)
    end
end

function tick_phage(phage, model)
    model[phage].time_in_state += 1
end

function infect(phage, cell, model)
    model[phage].state = :in_host
    model[phage].time_in_state = 0

    push!(model[cell].phages_inside, phage)
end

function burst(phage, cell, model)
    environment = model.properties.environment
    if environment === :well_mixed
        nearby = positions(model, :random)
    elseif environment === :spatially_structured
        r = 1
    elseif environment === :semi_solid
        r = 3
    end

    if environment === :spatially_structured || environment === :semi_solid
        nearby = collect(nearby_positions(model[phage].pos, model, r))
    end

    selected = rand(model.rng, nearby, model.properties.burst_size)
    for pos ∈ selected
        roll = rand(model.rng)
        if roll < 0.45
            kind = :temperate
        elseif roll < 0.85
            kind = :virulent
        else
            kind = :deficient
        end
        add_agent!(pos, Organism, model, roll < 0.5 ? :a : :b, Vector{Int}(), kind, :free, 0, :phage)
    end

    genocide!(model, model[cell].phages_inside)
    kill_agent!(cell, model)
end

function diffuse(id, model)
    agent = model[id]

    r = 3
    nearby = collect(nearby_positions(agent.pos, model, r))

    if agent.type === :phage
        selected = rand(model.rng, nearby)
        move_agent!(agent, selected, model)
    elseif agent.type === :bacterium
        move_bacterium_single!(agent, model, nearby)
    end
end

function grow(cell, p_grow, model)
    function has_no_bacteria(pos)
        ids = ids_in_position(pos, model)
        return isempty(ids) || !any(id -> model[id].type === :bacterium, ids)
    end

    !(rand(model.rng) < p_grow) && return

    environment = model.properties.environment
    if environment === :well_mixed
        selected = positions(model, :random)
    elseif environment === :spatially_structured
        r = 1
    elseif environment === :semi_solid
        r = 3
    end

    if environment === :spatially_structured || environment === :semi_solid
        selected = collect(nearby_positions(model[cell].pos, model, r))
    end

    filter!(has_no_bacteria, selected)
    isempty(selected) && return

    target = rand(model.rng, selected)
    add_agent!(target, Organism, model,
        rand(model.rng) < 0.5 ? :a : :b, Vector{Int}(), :nothing, :nothing, -1, :bacterium)
end

function bacteria_death_inherent(bacteria, model)
    for id ∈ bacteria
        if rand(model.rng) < p_death(model.a, model.b, model.m)
            genocide!(model, model[id].phages_inside)
            kill_agent!(id, model)
        end
    end
end

function bacteria_death_lysis(bacteria, model) # TODO: maybe only call this on infected bacteria (check the other (old) methods too)
    for cell ∈ bacteria
        phages_inside = model[cell].phages_inside
        if !isempty(phages_inside)
            for phage ∈ phages_inside
                tick_phage(phage, model)
                agent = model[phage]
                if (agent.time_in_state ≥ model.latent_period) && (agent.kind === :virulent || agent.kind === :induced_temperate)
                    if rand(model.rng) < model.properties.p_burst
                        burst(phage, cell, model)
                        break
                    end
                end
            end
        end
    end
end