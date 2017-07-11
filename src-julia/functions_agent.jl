###########
#Functions updating Agents
##########
function agent_age_factor(x::Int)
  f = sin(pi*x/AGE_MAX)
  return f
end

function agent_add_year!(a::Agent, t::Int)
  a._age += t
  agent_update!(a)
end
function agent_add_year!(p::Population, a::Agent, t::Int)
  a._age += t
  #agent_update!(p, a)

##
  if p.agents_random_replaced==true
    agent_update!(p, a)
  elseif p.agents_random_replaced==false
    agent_update_similar!(p, a)
  end
##
end

function agent_update!(a::Agent)
  a.age = a._age%AGE_MAX
  if a.age==0
    #if dead, replace with new agent with age 1
    a._age+=1
    a.age = a._age%AGE_MAX

    a._alpha = distribution(1)[1]
    a._beta = distribution(1)[1]
  end
  age_factor = agent_age_factor(a.age)
  a.alpha = a._alpha*age_factor
  a.beta = a._beta*age_factor
  a.dom_score = a.alpha+a.beta
end
function agent_update!(p::Population, a::Agent)
  a.age = a._age%AGE_MAX
  if a.age==0
    #if dead, replace with new agent with age 1
    a._age+=1
    a.age = a._age%AGE_MAX

    a._alpha = p._distribution_agents_params(1)[1]
    a._beta = p._distribution_agents_params(1)[1]
  end
  age_factor = agent_age_factor(a.age)
  a.alpha = a._alpha*age_factor
  a.beta = a._beta*age_factor
  a.dom_score = a.alpha+a.beta
end

function agent_update_similar!(a::Agent)
  a.age = a._age%AGE_MAX
  if a.age==0
    #if dead, replace with new agent with age 1
    a._age+=1
    a.age = a._age%AGE_MAX
  end
  age_factor = agent_age_factor(a.age)
  a.alpha = a._alpha*age_factor
  a.beta = a._beta*age_factor
  a.dom_score = a.alpha+a.beta
end
function agent_update_similar!(p::Population, a::Agent)
  a.age = a._age%AGE_MAX
  if a.age==0
    #if dead, replace with new agent with age 1
    a._age+=1
    a.age = a._age%AGE_MAX
  end
  age_factor = agent_age_factor(a.age)
  a.alpha = a._alpha*age_factor
  a.beta = a._beta*age_factor
  a.dom_score = a.alpha+a.beta
end

##################
################## Working on it
#################

##########
## Strategy: make function to compute indiv current best scores

function agent_get_current_dom_scores(p::Population, a::Agent)
  current_ds = population_get_all_groups_dom_scores(p)
  indices,groups_with_agent = population_get_groups_and_indices_containing_agent(p, a)
  agent_ds = current_ds[indices]
  return agent_ds
end
function agent_get_current_dom_scores(p::Population, a::Int)
  current_ds = population_get_all_groups_dom_scores(p)
  indices,groups_with_agent = population_get_groups_and_indices_containing_agent(p, a)
  agent_ds = current_ds[indices]
  return agent_ds
end

function agent_get_current_max_score(p::Population, a::Int)
  current_ds = population_get_all_groups_dom_scores(p)
  indices,groups_with_agent = population_get_groups_and_indices_containing_agent(p, a)
  agent_ds = current_ds[indices]
  return maximum(agent_ds)
end
function agent_get_current_max_score(p::Population, a::Agent)
  current_ds = population_get_all_groups_dom_scores(p)
  indices,groups_with_agent = population_get_groups_and_indices_containing_agent(p, a)
  agent_ds = current_ds[indices]
  return maximum(agent_ds)
end

##@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

function agent_random_coalitions(p::Population, a::Agent)
  #current dom scores (and max dom score) of the agent
  current_dom_scores = agent_get_current_dom_scores(p, a)
  current_max_dom_scores = maximum(current_dom_scores)

  #list of all the other agents
  other_agents = filter(ag->ag!=a, p.agents_list)
  nb_other_agents = length(other_agents)

  #random order of agents to make group with
  rdm_order = shuffle(collect(1:nb_other_agents))
  growing_groups = [push!(rdm_order[1:i], a.id) for i in 1:nb_other_agents]

  ds_list = [0. for i in 1:nb_other_agents]
  for i in 1:nb_other_agents
    agents = get_agents(p.agents_list, growing_groups[i])
    ds = get_dom_score_agents_list(agents, p.g_max, p.g_kd, p.g_h) ##this function changed to add laziness
    ds_list[i] = ds
  end
  return ds_list
end





function agent_random_coalition_with_groups(p::Population, a::Agent)
  #random order of groups indices
  rdm_order = shuffle(collect(1:length(p.groups_list)))

  #current list of agents per group in the population and respective dom scores
  current_agents_list = [Set(filter(a->a.id in g.agents_id_list, p.agents_list)) for g in p.groups_list]
  current_ds = population_get_all_groups_dom_scores(p)

  #growing list of agents
  new_agents_list = []
  new_ds = []
  growing_list = Set([a])
  for (i,agents) in enumerate(current_agents_list)
    union!(growing_list, Set(agents))
    push!(new_agents_list, growing_list)
    ds = get_dom_score_agents_list(growing_list, p) ##this function changed to add laziness
    push!(new_ds, ds)
  end
  return new_agents_list, new_ds
end

function agent_random_coalition_with_groups2(p::Population, a::Agent)
  #random order of groups indices
  rdm_order = shuffle(collect(1:length(p.groups_list)))

  #current list of agents per group in the population and respective dom scores
  current_agents_id_list = [g.agents_id_list for g in p.groups_list]
  current_ds = population_get_all_groups_dom_scores(p)

  #growing list of agents
  new_agents_list = copy(current_agents_id_list)
  push!(new_agents_list[1], a.id)
  new_ds = []
  growing_list = Set([a.id])
  for (i,agents) in enumerate(current_agents_list)
    union!(growing_list, Set(agents))
    push!(new_agents_list, growing_list)
    ds = get_dom_score_agents_list(growing_list, p) ##this function changed to add laziness
    push!(new_ds, ds)
  end
  return new_agents_list, new_ds
end
