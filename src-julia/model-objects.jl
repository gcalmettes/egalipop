using Distributions

###########
#Individuals
###########

type Agent
  id::Int
  group::Int
  _age::Int
  _alpha::Float64
  _beta::Float64
  age::Int
  alpha::Float64
  beta::Float64
  dom_score::Float64

  function Agent(idx::Int)
    age = rand(1:AGE_MAX-1)
    alpha = distribution(1)[1]
    beta = distribution(1)[1]
    alpha_age = alpha*agent_age_factor(age)
    beta_age = beta*agent_age_factor(age)
    new(idx, idx, age, alpha, beta, age%AGE_MAX, alpha_age, beta_age, alpha_age+beta_age)
  end
end


###########
#Population
###########

#abstract type so that circular definition of the Group type inside Population work
abstract AbstractGroup

type Population
  size::Int
  g_max::Float64
  g_h::Float64
  g_kd::Float64
  agents_list::Vector{Agent}
  groups_list::Vector{AbstractGroup} #::Vector{Group} doesn't work as circular definition is not yet handled by JULIA
  #groups_list::Vector{Group}
  _agents_count::Int #use to attribute new agents id
  _groups_count::Int #use to attribute new group id
  _distribution_agents_params::Function
  _distribution_groups_laziness::Distribution

  function Population(size::Int=15 ; gmax::Float64=G_MAX, gh::Float64=G_H, gkd::Float64=G_KD, distribution_agents_params::Function=dist_poisson, distribution_groups_laziness::Distribution=Uniform(0.5, 1), in_groups::Bool=false)
    p = new(size, gmax, gh, gkd, [Agent(i) for i in 1:size], Group[], size, 0, distribution_agents_params, distribution_groups_laziness)
    if in_groups==true
        for agent in p.agents_list
            population_add_groups!(p, 1, [agent])
        end
    end
    return p
  end
end


###########
#Groups
###########

type Group <: AbstractGroup
  id::Int
  population_linked_to::Population
  agents_id_list::Set{Int}
  laziness_value::Float64

  function Group(id::Int, p::Population)
    grp = new(id, p, Set(Int64[]), 1.)
    return grp
  end
  function Group(id::Int, p::Population, agents::Vector{Agent})
    agent_idx = get_agent_ids(agents)
    if length(agent_idx)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(id, p, Set(agent_idx), lambda)
    return grp
  end
  function Group(id::Int, p::Population, agents_id::Vector{Int})
    if length(agents_id)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(id, p, Set(agents_id), lambda)
    return grp
  end
  function Group(id::Int, p::Population, agents_id::Set{Int})
    if length(agents_id)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(id, p, agents_id, lambda)
    return grp
  end
  function Group(p::Population)
    grp = new(p._groups_count+1, p, Set(Int64[]), 1.)
    return grp
  end
  function Group(p::Population, agents::Vector{Agent})
    agent_idx = get_agent_ids(agents)
    if length(agent_idx)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(p._groups_count+1, p, Set(agent_idx), lambda)
    return grp
  end
  function Group(p::Population, agents_id::Vector{Int})
    if length(agents_id)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(p._groups_count+1, p, Set(agent_id), lambda)
    return grp
  end
  function Group(p::Population, agents_id::Set{Int})
    if length(agents_id)==1
      lambda = 1.
    else
      lambda = rand(p._distribution_groups_laziness)
    end
    grp = new(p._groups_count+1, p, agent_id, lambda)
    return grp
  end
  function Group(id::Int, p::Population, agents::Vector{Agent}, lambda::Float64)
    agent_idx = get_agent_ids(agents)
    grp = new(id, p, Set(agent_idx), lambda)
    return grp
  end
  function Group(id::Int, p::Population, agents_id::Vector{Int}, lambda::Float64)
    grp = new(id, p, Set(agents_id), lambda)
    return grp
  end
  function Group(id::Int, p::Population, agents_id::Set{Int}, lambda::Float64)
    grp = new(id, p, agents_id, lambda)
    return grp
  end
  function Group(p::Population, agents::Vector{Agent}, lambda::Float64)
    agent_idx = get_agent_ids(agents)
    grp = new(p._groups_count+1, p, Set(agent_idx), lambda)
    return grp
  end
  function Group(p::Population, agents_id::Vector{Int}, lambda::Float64)
    grp = new(p._groups_count+1, p, Set(agent_id), lambda)
    return grp
  end
  function Group(p::Population, agents_id::Set{Int}, lambda::Float64)
    grp = new(p._groups_count+1, p, agent_id, lambda)
    return grp
  end
end
