@agent Organism GridAgent{2} begin
    species::Symbol
    phages_inside::Vector{Int}
    kind::Symbol
    state::Symbol
    time_in_state::Int
    type::Symbol
end

Base.@kwdef mutable struct Parameters
    bacteria_count::Int = 265
    phages_count::Int = 148
    environment::Symbol = :semi_solid
    diffuse::Bool = true
    a::Float64 = 0.0
    b::Float64 = 0.08
    m::Float64 = 27.0
    α::Float64 = 8.0
    κ::Float64 = 0.05
    moi_proxy_radius::Int = 1
    infection_distance::Int = 1
    latent_period::Int = 5
    burst_size::Int = 4
    carrying_capacity::Int = 3 * bacteria_count
    growth_rate::Float64 = 0.6
    decay_factor::Float64 = 0.02
    p_burst = 0.9
end