# Bacteria-Phage Population Coevolution Dynamics ABM


## Description

This is an Agent-Based model of the coevolution dynamics between a bacteriaphage and bacteria population. It is written in [Julia](https://julialang.org/) and uses [`Agents.jl`](https://juliadynamics.github.io/Agents.jl/stable/) from [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/) for the built-in ABM logic. The model, as is usually the case when dealing with ABMs, consists of several simple equations, with deterministic and stochastic components. The goal is to explore a relatively new insight in the field: how the environment structure alters the coevolution dynamics. It draws inspiration by [^1], as published in [_Nature_](https://www.nature.com/articles/s41598-019-39773-3). The model is completely interactive, with a GUI that allows the user to modify the simulation parameters on the fly, as well as quickly visualize the dynamics with real-time plots. It is left simple yet complete, to encourage the user to further extend upon both the model design, but also the framework functionality (e.g. add more plots, automatically extract data, automatically run parametrized simulations...). The environment structure affects how easy it is for the cells and phage particles to change position in space. This mainly affects the diffusion (or lack thereof) process for both the cells and phage particles that are free in the environment, as well as the newly introduced phage particles post-cell burst. It also affects the mechanics that govern how the new cells are introduced spatially, etc.

[^1]: Sousa, J. A., & Rocha, E. P. (2019). Environmental structure drives resistance to phages and antibiotics during phage therapy and to invading lysogens during colonisation. Scientific reports, 9(1), 1-13.

## Table of Contents

<details>
<summary>Click to expand</summary>

- [Bacteria-Phage Population Coevolution Dynamics ABM](#bacteria-phage-population-coevolution-dynamics-abm)
  - [Description](#description)
  - [Table of Contents](#table-of-contents)
  - [Motivation](#motivation)
  - [Model Structure](#model-structure)
    - [Bacterial Death](#bacterial-death)
    - [Phage Decay](#phage-decay)
    - [Lysogeny](#lysogeny)
    - [Adsorption](#adsorption)
    - [Growth](#growth)
    - [](#)
  - [Implementation](#implementation)
  - [Gallery](#gallery)

</details>

## Motivation

<details>
<summary>Click to expand</summary>

Different kinds of what can be classified as a microbial organism can be found in every natural, or natural-alike, environment; the human body not being an exception. The populations of microbial organisms are key factors in how the environments they inhabit are temporally shaped. Not only do said organisms contribute in forming the greater enclosing ecosystem they come to be a part of, but they can also significantly affect it in various ways, the most important of which being forcing the host to a “reaction-chain” of constant adaptation, thus driving evolution. A long-standing study focus, rising in popularity, is in regard to the interactions between bacteria and bacteriophages (or simply, phages). The phages predate on the bacteria, thereby having a regulatory role in the bacteria growth and overall population dynamics, while also forming the “rules” by which they adapt. Phages are both incredible predators, as well as the most abundant entities in nature. Pathogenic bacteria showcase an increasingly effective resistance to antibiotics. Furthermore, the development of new antibiotics has slowed down considerably, not being able to keep pace with the surge in the appearance of more potent, unaffected by current treatment, bacteria. Thus, there has been rekindled interest in utilizing phages in a controlled setting to provide an alternative of or complement to an antibiotic treatment. As the progress in this field can be considered at a nascent stage, modeling the bacteria-phage interactions is of paramount importance towards exploiting them to treat disease.

In brief, phages are infectious acellular entities that depend on the existence of bacterial cells to proliferate. In typical predator-prey fashion, the phage population growth depends on the respective bacteria population growth. The host bacteria evolve mechanisms to resist the phages, while, at the same time, the phage particles evolve new strategies to manage to infect them, in an asynchronous manner. The huge variety of existing bacteria and phage strains, the vast number of discrete states, the non-linearity of the pharmacodynamics and pharmacokinetics in play, coupled with the inherent stochasticity that characterizes evolution, among others, render the bacteria-phage coevolution system encompassing complex dynamics.

To elucidate the mechanisms underlying phage-bacteria interactions, a variety of experimental in vitro and in vivo approaches have been developed, on simplified and more complex environments; natural along with simulated ones. Mathematical modeling contributes in propelling all related research forward, as well as helps combat some inevitably arising technicalities. For example, there is limited knowledge of the antagonistic coevolution in nature, it is difficult to maintain target and control cultures in parallel to carry out experiments, and only few bacteria are even amenable to being cultured in a laboratory. While, through modeling, we can arrive to analytical solutions derived from using well-established techniques to explore parameter and solution space, the models often don’t scale well, failing to address spatial heterogeneity (very important since the dynamics are heavily affected by environment structure) and key processes that are stochastic in nature. The models occasionally fail to pinpoint the effects of individual mechanisms and highly-specialized cases, have limited resolution in tracking temporal dynamics and may end up being intractable, as the system under study increases in complexity. These are some of the reasons why the more popular modeling approaches (using delay differential equations, cellular automata, MCMC, etc.) can be unable to reproduce the observational data, and/or not be transferable to realistic scenarios.

The above can pave the way for an agent-based modeling approach to incorporate different mechanisms at the level of the individual, supporting local interactions. It can constitute a way to include low-level biological detail, while having the system-level dynamics that emerge from the local, independent interactions and decisions remain intact and easily accessible.

</details>

## Model Structure

The model draws from many submodels that constitute mainstream approaches in modeling the coevolution dynamics between these two populations. The model space selected is a 2D Grid Space with Chebyshev metric.

### Bacterial Death

```julia
p_death(a, b, m) = a + ((1 - a) / (1 + exp(-b * (-m))))
```

$$ DeathProbability = A + \frac{1-A}{1+e^{-B(-M)}} $$

Where `A`, `B` and `M` are three parameters to control the shape of the death curve.

### Phage Decay

```julia
p_phage_decay(decay_factor, time_in_state) = 1 - (1 * exp(-decay_factor * time_in_state))
```

$$ PhageDecayProbability = 1 - (1 \cdot e^{-DecayFactor \cdot TimeOutsideHost}) $$

Where `DecayFactor` controls the environmental degradation rate, while `TimeOutsideHost` is the number of generation the particiular phage particle has spent free in the environment (i.e. not "inside" a cell host).

### Lysogeny

```julia
function p_lysis(phage, model)
    nearby_phages = nearby_t(:phage, phage, model)
    nearby_phages = length(nearby_phages)
    return 1 / (1 + model.properties.α * exp(-nearby_phages + model.properties.κ))
end
```

$$ LysogenyProbability = \frac{1}{1 + \alpha \cdot e^{-SurroundingPhages + \kappa}} $$

Where `SurroundingPhages` is the number of phage particles within Moore distance dictated by the infection distance (see below), while the parameters `α` and `κ` control the lysogeny curve.

### Adsorption

```julia
function p_adsorption(cell, model)
    nearby_phages = nearby_t(:phage, cell, model)
    p_host = model[cell].species === :a ? 0.9 : 0.1
    return p_host / (1 + exp(-length(nearby_phages)))
end
```

$$ AdsorptionProbability = \frac{HostAdsorptionProbability}{1 + e^{SurroundingPhages}} $$

### Growth

```julia
function p_grow(model)
    cells = by_single_type(:bacterium)(model)
    filter!(id -> isempty(model[id].phages_inside), cells)
    cell_count = length(cells)

    properties = model.properties
    ΔN = properties.growth_rate * cell_count * (1 - cell_count / properties.carrying_capacity)
    return ΔN / cell_count
end
```

$$ DivisionProbability = \frac{\Delta N}{CellCount} $$

Where

$$ \Delta N  = GrowthRate \cdot CellCount (1 - \frac{CellCount}{CarryingCapacity})$$

A simple implementation of the logistic equation as introduced by Verhulst. Governs whether the cell will divide, adding a new cell in the population.

### 

## Implementation

In a Monte Carlo method fashion, the generated probabilities (e.g. for cell division, death, phage decay...) are used for repeated random sampling.

For each generation:
1. Cells die according to inherent causes
2. Infected cells die due to lysis, new phage particles are introduced  to the environment as the cellular structure bursts
3. For all phage particles that are free in the environment, `TimeOutsideHost` gets updated
4. Free phage particles infect cells according to the adsorption probability
5. Free phage particles decay due to environmental degradation
6. Uninfected cells grow (divide)
7. According to the environment, cells and phage particles diffuse

```julia
function complex_step!(model)
    bacteria_death_inherent(by_single_type(:bacterium)(model), model)
    bacteria_death_lysis(by_single_type(:bacterium)(model), model)

    phages = by_single_type(:phage)(model)
    if isempty(phages)
        (model.properties.phages_count == 0) && return

        model.properties.bacteria_count = length(by_single_type(:bacterium)(model))
        model.properties.phages_count = 0
        return
    end

    filter!(id -> model[id].state === :free, phages)
    for phage ∈ phages
        free_phage_step(phage, model)
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
        cells = by_single_type(:bacterium)(model)
        phages = by_single_type(:phage)(model)
        filter!(id -> model[id].state === :free, phages)
        if model.properties.environment === :semi_solid
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
```

Where

```julia
function free_phage_step(phage, model)
    agent = model[phage]

    nearby_cells = nearby_t(:bacterium, phage, model)
    isempty(nearby_cells) && return

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
```

## Gallery

![Initial stage](docs/assets/abm1.png)

---

![Progression 1](docs/assets/abm2.png)

---

![Progression 2](docs/assets/abm3.png)

---

![Progression 3](docs/assets/abm4.png)

---

![Progression 4](docs/assets/abm5.png)
