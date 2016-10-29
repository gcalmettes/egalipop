##########
#Functions related to Group object
##########
function group_get_dom_score(g::Group)
  agents_list = filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list)
  ds = get_dom_score_agents_list(agents_list,
                                 g.population_linked_to.g_max,
                                 g.population_linked_to.g_kd,
                                 g.population_linked_to.g_h,
                                 g.laziness_value)
  return ds
end
function group_get_dom_score(p::Population, g_id::Int)
  g = population_get_group(p, g_id)
  agents_list = filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list)
  ds = get_dom_score_agents_list(agents_list,
                                 g.population_linked_to.g_max,
                                 g.population_linked_to.g_kd,
                                 g.population_linked_to.g_h,
                                 g.laziness_value)
  return ds
end

function group_get_size(g::Group)
  s = length(g.agents_id_list)
  return s
end
function group_get_size(p::Population, g_id::Int)
  g = population_get_group(p, g_id)
  s = length(g.agents_id_list)
  return s
end

function group_add_agents!(g::Group, a::Agent)
  push!(g.agents_id_list, a.id)
end
function group_add_agents!(g::Group, a::Int)
  push!(g.agents_id_list, a)
end
function group_add_agents!(g::Group, a::Vector{Agent})
  new_ids = Set([ag.id for ag in a])
  union!(g.agents_id_list, new_ids)
end
function group_add_agents!(g::Group, a::Vector{Int})
  union!(g.agents_id_list, Set(a))
end
function group_add_agents!(g::Group, a::Set{Int})
  union!(g.agents_id_list, a)
end
function group_add_agents!(p::Population, g_id::Int, a::Agent)
  g = population_get_group(p, g_id)
  push!(g.agents_id_list, a.id)
end
function group_add_agents!(p::Population, g_id::Int, a::Int)
  g = population_get_group(p, g_id)
  push!(g.agents_id_list, a)
end
function group_add_agents!(p::Population, g_id::Int, a::Vector{Agent})
  g = population_get_group(p, g_id)
  new_ids = Set([ag.id for ag in a])
  union!(g.agents_id_list, new_ids)
end
function group_add_agents!(p::Population, g_id::Int, a::Vector{Int})
  g = population_get_group(p, g_id)
  union!(g.agents_id_list, Set(a))
end
function group_add_agents!(p::Population, g_id::Int, a::Set{Int})
  g = population_get_group(p, g_id)
  union!(g.agents_id_list, a)
end

function group_remove_agents!(g::Group, a::Int)
  filter!(x->x!=a, g.agents_id_list)
end
function group_remove_agents!(g::Group, a::Agent)
  filter!(x->x!=a.id, g.agents_id_list)
end
function group_remove_agents!(g::Group, a::Vector{Int})
  for id in a
    filter!(x->x!=id, g.agents_id_list)
  end
end
function group_remove_agents!(g::Group, a::Vector{Agent})
  ag_ids = [ag.id for ag in a]
  for id in ag_ids
    filter!(x->x!=id, g.agents_id_list)
  end
end
function group_remove_agents!(p::Population, g_id::Int, a::Int)
  g = population_get_group(p, g_id)
  filter!(x->x!=a, g.agents_id_list)
end
function group_remove_agents!(p::Population, g_id::Int, a::Agent)
  g = population_get_group(p, g_id)
  filter!(x->x!=a.id, g.agents_id_list)
end
function group_remove_agents!(p::Population, g_id::Int, a::Vector{Int})
  g = population_get_group(p, g_id)
  for id in a
    filter!(x->x!=id, g.agents_id_list)
  end
end
function group_remove_agents!(p::Population, g_id::Int, a::Vector{Agent})
  g = population_get_group(p, g_id)
  ag_ids = [ag.id for ag in a]
  for id in ag_ids
    filter!(x->x!=id, g.agents_id_list)
  end
end

function group_dom_score_with_agent(g::Group, a::Agent)
  ###Return dom score of group g if agent a was part of the group
  agents_list = Set(filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list))
  push!(agents_list, a)
  new_ds = get_dom_score_agents_list(agents_list, g.population_linked_to, g.laziness_value)
  return new_ds
end
function group_dom_score_with_agent(p::Population, g::Int, a::Int)
  ###Return dom score of group g if agent a was part of the group
  g = population_get_group(p, g)
  a = population_get_agent(p, a)
  agents_list = Set(filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list))
  push!(agents_list, a)
  new_ds = get_dom_score_agents_list(agents_list, g.population_linked_to, g.laziness_value)
  return new_ds
end
function group_dom_score_with_agent(g::Group, a::Agent, lambda::Float64)
  ###Return dom score of group g if agent a was part of the group
  agents_list = Set(filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list))
  push!(agents_list, a)
  new_ds = get_dom_score_agents_list(agents_list, g.population_linked_to, lambda)
  return new_ds
end
function group_dom_score_with_agent(p::Population, g::Int, a::Int, lambda::Float64)
  ###Return dom score of group g if agent a was part of the group
  g = population_get_group(p, g)
  a = population_get_agent(p, a)
  agents_list = Set(filter(a->a.id in g.agents_id_list, g.population_linked_to.agents_list))
  push!(agents_list, a)
  new_ds = get_dom_score_agents_list(agents_list, g.population_linked_to, lambda)
  return new_ds
end
