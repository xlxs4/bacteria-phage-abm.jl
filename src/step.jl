function complex_step!(model)
    bacteria_death_inherent(by_single_type(:bacterium)(model), model)
    bacteria_death_lysis(by_single_type(:bacterium)(model), model)

    phages = by_single_type(:phage)(model)
    if isempty(phages)
        (model.properties.phages_count == 0) && return nothing

        model.properties.bacteria_count = length(by_single_type(:bacterium)(model))
        model.properties.phages_count = 0
        return nothing
    end

    filter!(id -> model[id].state === :free, phages)
    for phage ∈ phages
        agent = model[phage]

        nearby_cells = nearby_t(:bacterium, phage, model)
        isempty(nearby_cells) && continue

        target_cell = rand(model.rng, nearby_cells)
        if rand(model.rng) < p_adsorption(target_cell, model)
            kind = agent.kind
            if kind === :temperate
                if rand(model.rng) < p_lysis(phage, model)
                    agent.kind = :induced_temperate
                end
            end
            if kind === :virulent || kind === :induced_temperate
                infect(phage, target_cell, model)
            end
        end
    end

    filter!(id -> model[id].state === :free, phages)
    for phage ∈ phages
        tick_phage(phage, model)
        phage_decay(phage, model)
    end

    cells = by_single_type(:bacterium)(model)
    filter!(id -> isempty(model[id].phages_inside), cells)
    p = p_grow(model)
    for cell ∈ cells
        grow(cell, p, model)
    end

    if model.properties.diffuse
        environment = model.properties.environment
        cells = by_single_type(:bacterium)(model)
        phages = by_single_type(:phage)(model)
        filter!(id -> model[id].state === :free, phages)
        if environment === :semi_solid
            for cell ∈ cells
                diffuse(cell, model)
            end
            for phage ∈ phages
                diffuse(phage, model)
            end
        end
    end


    model.properties.bacteria_count = length(by_single_type(:bacterium)(model))
    model.properties.phages_count = length(by_single_type(:phage)(model))
end