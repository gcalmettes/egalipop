using Iterators

###########
#Functions related to Population object
##########
function population_get_alphas(p::Population)
  alphas = [a.alpha for a in p.agents_list]
end

function population_get_betas(p::Population)
  betas = [a.beta for a in p.agents_list]
end

function population_get_agent(p::Population, id::Int)
  agent = filter(a->a.id==id, p.agents_list)
  if length(agent)==0
    print("Agent ", id, " doesn't exist")
    return []
  else
    return agent[1]
  end
end

function population_get_group(p::Population, id::Int)
  group = filter(g->g.id==id, p.groups_list)
  return group[1]
end

function population_get_group_dom_score(p::Population, id::Int)
  g = population_get_group(p, id)
  agents_list = filter(a->a.id in g.agents_id_list, p.agents_list)
  ds = get_dom_score_agents_list(agents_list, p.g_max, p.g_kd, p.g_h, g.laziness_value)
  return ds
end

function population_get_all_groups_dom_scores_with_agent(p::Population, a::Agent)
  ###Return dom score of all the groups if agent "a" was
  ###part of each group
  agents_per_group = [Set(filter(a->a.id in g.agents_id_list, p.agents_list)) for g in p.groups_list]
  groups_with_agent = population_get_groups_containing_agent(p, a)
  laziness_per_group = [g.laziness_value for g in groups_with_agent]
  nb_groups = length(agents_per_group)
  new_ds_per_group = Array(Float64, nb_groups)
  for i = 1:nb_groups
    push!(agents_per_group[i], a)
    new_ds = get_dom_score_agents_list(agents_per_group[i], p, laziness_per_group[i])
    new_ds_per_group[i]=new_ds
  end
  return new_ds_per_group
end
function population_get_all_groups_dom_scores_with_agent(p::Population, a::Int)
  ###Return dom score of all the groups if agent "a" was
  ###part of each group
  a = population_get_agent(p, a)
  agents_per_group = [Set(filter(a->a.id in g.agents_id_list, p.agents_list)) for g in p.groups_list]
  groups_with_agent = population_get_groups_containing_agent(p, a)
  laziness_per_group = [g.laziness_value for g in groups_with_agent]
  nb_groups = length(agents_per_group)
  new_ds_per_group = Array(Float64, nb_groups)
  for i = 1:nb_groups
    push!(agents_per_group[i], a)
    new_ds = get_dom_score_agents_list(agents_per_group[i], p, laziness_per_group[i])
    new_ds_per_group[i]=new_ds
  end
  return new_ds_per_group
end

function population_get_group_size(p::Population, id::Int)
  g = population_get_group(p, id)
  size = length(g.agents_id_list)
  return size
end

function population_get_all_groups_dom_scores(p::Population)
  groups_count = length(p.groups_list)
  ds = Array(Float64, groups_count)
  for i in 1:groups_count
    g = p.groups_list[i]
    agents_list = filter(a->a.id in g.agents_id_list, p.agents_list)
    g_ds = get_dom_score_agents_list(agents_list, p.g_max, p.g_kd, p.g_h, g.laziness_value)
    ds[i] = g_ds
  end
  return ds
end

function population_get_all_agents_dom_scores(p::Population)
  size = p.size
  ds = Array(Float64, size)
  for i in 1:size
    ds[i] = agent_get_current_max_score(p, i)
  end
  return ds
end

function population_get_all_groups_sizes(p::Population)
  s = [length(g.agents_id_list) for g in p.groups_list]
  return s
end

function population_get_non_empty_groups(p::Population)
  sizes = population_get_all_groups_sizes(p)
  valid_groups = []
  for i in 1:length(p.groups_list)
    if sizes[i]>0
      push!(valid_groups, p.groups_list[i])
    end
  end
  return valid_groups
end

function population_add_agents!(p::Population, n::Int)
  start_id = p._agents_count
  new_agents = [Agent(i) for i in start_id+1:start_id+n]
  append!(p.agents_list, new_agents)
  population_update_agents_count!(p, n)
  p.size += n
end

function population_add_groups!(p::Population, n::Int)
  new_groups = [Group(p._groups_count+i, p) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Agent})
  new_groups = [Group(p._groups_count+i, p, a) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Int})
  new_groups = [Group(p._groups_count+i, p, a) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, lambda::Float64)
  new_groups = [Group(p._groups_count+i, p, lambda) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Agent}, lambda::Float64)
  new_groups = [Group(p._groups_count+i, p, a, lambda) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Int}, lambda::Float64)
  new_groups = [Group(p._groups_count+i, p, a, lambda) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, lambda::Vector{Float64})
  new_groups = [Group(p._groups_count+i, p, lambda[i]) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Agent}, lambda::Vector{Float64})
  new_groups = [Group(p._groups_count+i, p, a, lambda[i]) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end
function population_add_groups!(p::Population, n::Int, a::Vector{Int}, lambda::Vector{Float64})
  new_groups = [Group(p._groups_count+i, p, a, lambda[i]) for i in 1:n]
  append!(p.groups_list, new_groups)
  population_update_groups_count!(p, n)
end

function population_update_groups_count!(p::Population)
  p._groups_count += 1
end
function population_update_groups_count!(p::Population, i::Int)
  p._groups_count += i
end

function population_update_agents_count!(p::Population)
  p._agents_count += 1
end
function population_update_agents_count!(p::Population, i::Int)
  p._agents_count += i
end

function population_generate_random_groups!(p::Population; nboot::Int=1000)
    for i in 1:nboot
        rdm_splits = random_split(p.agents_list)
        for alist in rdm_splits
            population_add_groups!(p, 1, alist)
        end
    end
end

function population_get_groups_containing_agent(p::Population, a::Agent)
    groups = []
    for g in p.groups_list
        if a.id in g.agents_id_list
            push!(grps, g)
        end
    end
  return groups
end
function population_get_groups_containing_agent(p::Population, a::Int)
    groups = []
    for g in p.groups_list
        if a in g.agents_id_list
            push!(grps, g)
        end
    end
  return groups
end

function population_get_groups_and_indices_containing_agent(p::Population, a::Agent)
  indices = []
  groups = []
  for (i,g) in enumerate(p.groups_list)
    if a.id in g.agents_id_list
      push!(indices, i)
      push!(groups, g)
    end
  end
  return indices,groups
end
function population_get_groups_and_indices_containing_agent(p::Population, a::Int)
  indices = []
  groups = []
  for (i,g) in enumerate(p.groups_list)
    if a in g.agents_id_list
      push!(indices, i)
      push!(groups, g)
    end
  end
  return indices, groups
end

function population_remove_agent!(p::Population, a::Int)
  agent = population_get_agent(p, a)
  filter!(x->x!=agent, p.agents_list)
  groups_to_purge = population_get_groups_containing_agent(p, a)
  for group in groups_to_purge
    group_remove_agents!(group, a)
  end
end
function population_remove_agent!(p::Population, a::Agent)
  filter!(x->x!=a, p.agents_list)
  groups_to_purge = population_get_groups_containing_agent(p, a)
  for group in groups_to_purge
    group_remove_agents!(group, a)
  end
end

function population_apply_time_increment!(p::Population, t::Int)
  for agent in p.agents_list
    agent_add_year!(p, agent, t)
  end
end


############################
###########################

function population_get_all_agents_score_differences(p::Population)
    all_scores = population_get_all_agents_dom_scores(p)
    all_pairs = collect(combinations(all_scores, 2))

    n = length(all_pairs)
    differences = Array(Float64, n)
    for i = 1:n
        a,b = all_pairs[i]
        differences[i] = abs(a-b)
    end
    return differences
end


############################
###########################


function reached_equilibrium(a::Vector{Float64}, b::Vector{Float64})
  return a==b
end

#=
function population_random_coalitions!(p::Population)

  pattern_end = [0. for i in 1:2]
  pattern_start = [1. for i in 1:2]

  while !reached_equilibrium(pattern_end, pattern_start) #loop until agents are no more changing groups
    pattern_start = population_get_all_groups_dom_scores(p)
    ds = copy(pattern_start) #population_get_all_groups_dom_scores(p)
    #non_empty_indices = find(x->x!=0., ds)
    glist = p.groups_list#[non_empty_indices]
    #ds = ds[non_empty_indices]
    current_groups = [(Float64(score), Int64(grp.id), grp) for (score, grp) in zip(ds, glist)] #score/id/group
    sort!(current_groups, by=x->x[1]) #sort by dom score

    nb_groups = length(current_groups)
    groups_agents_and_ids = [(shuffle(filter(a->a.id in g.agents_id_list, p.agents_list)), id) for (s,id,g) in current_groups]

    #iterate over the sorted list of groups
    for i = 1:nb_groups
      #tic()
      (target_group,target_group_id) = groups_agents_and_ids[i]
      if length(target_group)==0 #if empty group, iterate over next group
        continue
      end
      others_groups_and_ids = deleteat!(copy(groups_agents_and_ids), i)

      #############
      ##get rid of all but 1 empty group
      ##prevent doing several time the same calculation while still
      ##ensuring that agent will leave group if the individual score is higher
      ##than ds of any group he could be part of (since solo group will be checked)
      ############
      unique_combination = []
      pos_to_keep = []
      for m in 1:nb_groups-1
        g2check = others_groups_and_ids[m]
        if g2check[1] in unique_combination
          continue
        else
          push!(pos_to_keep, m)
        end
        push!(unique_combination, g2check[1])
      end
      others_groups_and_ids = [others_groups_and_ids[y] for y in pos_to_keep]



      ###### Note:
      ###### Evaluate option of doing the shuffle of groups to make alliances with
      ###### at each new agent, instead of each new groups

      shuffle!(others_groups_and_ids) #shuffle order of other groups to make alliances with

      #add target group at the end so the full population will also be checked
      #in case one of the agents of the target group is in none of the others groups
      push!(others_groups_and_ids, (target_group,target_group_id))

      others_groups = [g for (g,id) in others_groups_and_ids]
      others_ids = [id for (g,id) in others_groups_and_ids]
      #push!(others_groups, target_group)
      #push!(others_ids, target_group_id)
      others_groups_agents = vcat(others_groups...) #agents from random group order in 1 array

      gsizes = [length(grp) for grp in others_groups]
      sum_gsizes = cumsum(gsizes)

      #list of unique agents to form coalition with and IDs of groups of those agents
      groups_to_join = [Set(others_groups_agents[1:k]) for k in sum_gsizes]
      groups_ids_to_join = [Vector{Int}(others_ids[1:k]) for k in 1:length(gsizes)]
      #println(groups_ids_to_join)
      #each agent of target group tries to defect to others groups or groups of groups
      for agent in target_group
        #current_dom_score of agent
        current_dom_scores = agent_get_current_dom_scores(p, agent)
        current_max_dom_score = maximum(current_dom_scores)

        #scores of the other groups with the agent added
        potential_groups = [push!(copy(grp), agent) for grp in groups_to_join]
        potential_groups_ds = [get_dom_score_agents_list(list, p) for list in potential_groups]

        #check if better than current max score
        improved_scores = compare_float_vector(current_max_dom_score, potential_groups_ds) #this is a boolean array
        #improved_groups = potential_groups[improved_scores]

        nb_better_grp = sum(improved_scores)

        if nb_better_grp>0
          ##which of those groups provide improvement for all the members?
          groups_to_check = potential_groups[improved_scores]
          ds_to_check = potential_groups_ds[improved_scores]
          ids_to_check = groups_ids_to_join[improved_scores]
          better_groups = falses(nb_better_grp)
          #println(nb_better_grp, " ", length(groups_to_check))
          for j = 1:nb_better_grp
            #println("... ", j, "/", nb_better_grp)
            better_groups[j] = check_if_valid_group(p, groups_to_check[j], ds_to_check[j])
          end

          if sum(better_groups)>0
            valid_groups = groups_to_check[better_groups]
            valid_ds = ds_to_check[better_groups]
            valid_ids = ids_to_check[better_groups]
            #join maximum scoring group

            new_group_idx = find(x->x==maximum(valid_ds), valid_ds)[end]#[end]
            #new_group_idx = find(x->x==rand(valid_ds), valid_ds)[1]#[end]

            new_group = valid_groups[new_group_idx]
            ids_to_merge_with = valid_ids[new_group_idx]
            agent_ids = [a.id for a in new_group]
            #println("agent ", agent.id, " from group ", target_group_id,
            #        " would benefit of making a group ", agent_ids,
            #        " with agents from groups ", ids_to_merge_with)
            #remove agent from original group ...


            check1 = [grp.agents_id_list for grp in p.groups_list]####

            group_remove_agents!(p, target_group_id, agent)

            check2 = [grp.agents_id_list for grp in p.groups_list]####

            if sum(population_get_all_groups_sizes(p))>15###
              println(check1)
              println("... agent ", agent.id, " removed from group ", target_group_id)
              println(check2)
            end###


            #println("... agent ", agent.id, " removed from group ", target_group_id)
            #add agents to the biggest groups of the ids to merge with
            grps = filter(g->g.id in ids_to_merge_with, p.groups_list)
            gsizes = [length(g.agents_id_list) for g in grps]
            final_grp_id = grps[find(x->x==maximum(gsizes), gsizes)[1]].id
            group_add_agents!(p, final_grp_id, [agent])
            #println("... agent ", agent.id, " added to group ", final_grp_id)
            for gp in ids_to_merge_with
              if gp==final_grp_id
                continue
              else
                g_start = population_get_group(p, gp)
                group_add_agents!(p, final_grp_id, g_start.agents_id_list)
                #println("... agents ", g_start.agents_id_list, " added to group ", final_grp_id)
                g_start.agents_id_list = Set(Int64[]) #clear agent_id_list
                #println("... agents list of group ", g_start.id, " cleared")
              end
            end
          end
        end
      end
      #print(i, ' ')
      #toc()
    end
    pattern_end = population_get_all_groups_dom_scores(p)
  end
end



function population_random_coalitions_probability!(p::Population, prob::Float64=0.5)

  pattern_end = [0. for i in 1:2]
  pattern_start = [1. for i in 1:2]

  while !reached_equilibrium(pattern_end, pattern_start) #loop until agents are no more changing groups
    pattern_start = population_get_all_groups_dom_scores(p)
    ds = copy(pattern_start) #population_get_all_groups_dom_scores(p)
    #non_empty_indices = find(x->x!=0., ds)
    glist = p.groups_list#[non_empty_indices]
    #ds = ds[non_empty_indices]
    current_groups = [(Float64(score), Int64(grp.id), grp) for (score, grp) in zip(ds, glist)] #score/id/group
    sort!(current_groups, by=x->x[1]) #sort by dom score

    nb_groups = length(current_groups)
    groups_agents_and_ids = [(shuffle(filter(a->a.id in g.agents_id_list, p.agents_list)), id) for (s,id,g) in current_groups]

    #iterate over the sorted list of groups
    for i = 1:nb_groups
      #tic()
      (target_group,target_group_id) = groups_agents_and_ids[i]
      if length(target_group)==0 #if empty group, iterate over next group
        continue
      end
      others_groups_and_ids = deleteat!(copy(groups_agents_and_ids), i)

      #############
      ##get rid of all but 1 empty group
      ##prevent doing several time the same calculation while still
      ##ensuring that agent will leave group if the individual score is higher
      ##than ds of any group he could be part of (since solo group will be checked)
      ############
      unique_combination = []
      pos_to_keep = []
      for m in 1:nb_groups-1
        g2check = others_groups_and_ids[m]
        if g2check[1] in unique_combination
          continue
        else
          push!(pos_to_keep, m)
        end
        push!(unique_combination, g2check[1])
      end
      others_groups_and_ids = [others_groups_and_ids[y] for y in pos_to_keep]



      ###### Note:
      ###### Evaluate option of doing the shuffle of groups to make alliances with
      ###### at each new agent, instead of each new groups

      shuffle!(others_groups_and_ids) #shuffle order of other groups to make alliances with

      #add target group at the end so the full population will also be checked
      #in case one of the agents of the target group is in none of the others groups
      push!(others_groups_and_ids, (target_group,target_group_id))

      others_groups = [g for (g,id) in others_groups_and_ids]
      others_ids = [id for (g,id) in others_groups_and_ids]
      #push!(others_groups, target_group)
      #push!(others_ids, target_group_id)
      others_groups_agents = vcat(others_groups...) #agents from random group order in 1 array

      gsizes = [length(grp) for grp in others_groups]
      sum_gsizes = cumsum(gsizes)

      #list of unique agents to form coalition with and IDs of groups of those agents
      groups_to_join = [Set(others_groups_agents[1:k]) for k in sum_gsizes]
      groups_ids_to_join = [Vector{Int}(others_ids[1:k]) for k in 1:length(gsizes)]
      #println(groups_ids_to_join)
      #each agent of target group tries to defect to others groups or groups of groups
      for agent in target_group
        #current_dom_score of agent
        current_dom_scores = agent_get_current_dom_scores(p, agent)
        current_max_dom_score = maximum(current_dom_scores)

        #scores of the other groups with the agent added
        potential_groups = [push!(copy(grp), agent) for grp in groups_to_join]
        potential_groups_ds = [get_dom_score_agents_list(list, p) for list in potential_groups]

        #check if better than current max score
        improved_scores = compare_float_vector(current_max_dom_score, potential_groups_ds) #this is a boolean array
        #improved_groups = potential_groups[improved_scores]

        nb_better_grp = sum(improved_scores)

        if nb_better_grp>0
          ##which of those groups provide improvement for all the members?
          groups_to_check = potential_groups[improved_scores]
          ds_to_check = potential_groups_ds[improved_scores]
          ids_to_check = groups_ids_to_join[improved_scores]
          better_groups = falses(nb_better_grp)
          #println(nb_better_grp, " ", length(groups_to_check))
          for j = 1:nb_better_grp
            #println("... ", j, "/", nb_better_grp)
            better_groups[j] = check_if_valid_group(p, groups_to_check[j], ds_to_check[j])
          end

          if sum(better_groups)>0
            #draw random number and if below probability of group formation do it
            rdm_draw = rand()

            if rdm_draw<prob
              valid_groups = groups_to_check[better_groups]
              valid_ds = ds_to_check[better_groups]
              valid_ids = ids_to_check[better_groups]
              #join maximum scoring group

              new_group_idx = find(x->x==maximum(valid_ds), valid_ds)[end]#[end]
              #new_group_idx = find(x->x==rand(valid_ds), valid_ds)[1]#[end]

              new_group = valid_groups[new_group_idx]
              ids_to_merge_with = valid_ids[new_group_idx]
              agent_ids = [a.id for a in new_group]
              #println("agent ", agent.id, " from group ", target_group_id,
              #        " would benefit of making a group ", agent_ids,
              #        " with agents from groups ", ids_to_merge_with)
              #remove agent from original group ...
              group_remove_agents!(p, target_group_id, agent)
              #println("... agent ", agent.id, " removed from group ", target_group_id)
              #add agents to the biggest groups of the ids to merge with
              grps = filter(g->g.id in ids_to_merge_with, p.groups_list)
              gsizes = [length(g.agents_id_list) for g in grps]
              final_grp_id = grps[find(x->x==maximum(gsizes), gsizes)[1]].id
              group_add_agents!(p, final_grp_id, [agent])
              #println("... agent ", agent.id, " added to group ", final_grp_id)
              for gp in ids_to_merge_with
                if gp==final_grp_id
                  continue
                else
                  g_start = population_get_group(p, gp)
                  group_add_agents!(p, final_grp_id, g_start.agents_id_list)
                  #println("... agents ", g_start.agents_id_list, " added to group ", final_grp_id)
                  g_start.agents_id_list = Set(Int64[]) #clear agent_id_list
                  #println("... agents list of group ", g_start.id, " cleared")
                end
              end
            else:
              continue
            end
          end
        end
      end
      #print(i, ' ')
      #toc()
    end
    pattern_end = population_get_all_groups_dom_scores(p)
  end
end

=#


function population_random_coalitions_freerider!(p::Population)

  pattern_end = [0. for i in 1:2]
  pattern_start = [1. for i in 1:2]
  while !reached_equilibrium(pattern_end, pattern_start) #loop until agents are no more changing groups
    pattern_start = population_get_all_groups_dom_scores(p)
    ds = copy(pattern_start) #population_get_all_groups_dom_scores(p)
    #non_empty_indices = find(x->x!=0., ds)
    glist = p.groups_list#[non_empty_indices]
    #ds = ds[non_empty_indices]
    current_groups = [(Float64(score), Int64(grp.id), grp) for (score, grp) in zip(ds, glist)] #score/id/group
    sort!(current_groups, by=x->x[1]) #sort by dom score

    nb_groups = length(current_groups)
    groups_agents_and_ids = [(shuffle(filter(a->a.id in g.agents_id_list, p.agents_list)), id) for (s,id,g) in current_groups]

    movement = 0

    #iterate over the sorted list of groups
    for i = 1:nb_groups
      #tic()
      (target_group,target_group_id) = groups_agents_and_ids[i]
      if length(target_group)==0 #if empty group, iterate over next group
        continue
      end
      others_groups_and_ids = deleteat!(copy(groups_agents_and_ids), i)

      #############
      ##get rid of all but 1 empty group
      ##prevent doing several time the same calculation while still
      ##ensuring that agent will leave group if the individual score is higher
      ##than ds of any group he could be part of (since solo group will be checked)
      ############
      unique_combination = []
      pos_to_keep = []
      for m in 1:nb_groups-1
        g2check = others_groups_and_ids[m]
        if g2check[1] in unique_combination
          continue
        else
          push!(pos_to_keep, m)
        end
        push!(unique_combination, g2check[1])
      end
      others_groups_and_ids = [others_groups_and_ids[y] for y in pos_to_keep]

      ###### Note:
      ###### Evaluate option of doing the shuffle of groups to make alliances with
      ###### at each new agent, instead of each new groups

      shuffle!(others_groups_and_ids) #shuffle order of other groups to make alliances with

      #add target group at the end so the full population will also be checked
      #in case one of the agents of the target group is in none of the others groups
      push!(others_groups_and_ids, (target_group,target_group_id))

      others_groups = [g for (g,id) in others_groups_and_ids]
      others_ids = [id for (g,id) in others_groups_and_ids]
      #push!(others_groups, target_group)
      #push!(others_ids, target_group_id)
      others_groups_agents = vcat(others_groups...) #agents from random group order in 1 array

      gsizes = [length(grp) for grp in others_groups]
      sum_gsizes = cumsum(gsizes)

      #list of unique agents to form coalition with and IDs of groups of those agents
      groups_to_join = [Set(others_groups_agents[1:k]) for k in sum_gsizes]
      groups_ids_to_join = [Vector{Int}(others_ids[1:k]) for k in 1:length(gsizes)]
      #println(groups_ids_to_join)
      #each agent of target group tries to defect to others groups or groups of groups
      for agent in target_group
        #current_dom_score of agent
        current_dom_scores = agent_get_current_dom_scores(p, agent)
        current_max_dom_score = maximum(current_dom_scores)

        #scores of the other groups with the agent added
        potential_groups = [push!(copy(grp), agent) for grp in groups_to_join]
        ###### LAZINESS FACTOR for potential groups
        potential_lambdas = rand(p._distribution_groups_laziness, length(potential_groups))
        for l in 1:length(potential_groups)
          if length(potential_groups[l])==1
            potential_lambdas[l]=1.
          end
        end
        #corresponding DS for potential groups
        potential_groups_ds = [get_dom_score_agents_list(list, p, lambda) for (list, lambda) in zip(potential_groups, potential_lambdas)]


        #check if better than current max score
        improved_scores = compare_float_vector(current_max_dom_score, potential_groups_ds) #this is a boolean array
        #improved_groups = potential_groups[improved_scores]

        nb_better_grp = sum(improved_scores)

        if nb_better_grp>0
          ##which of those groups provide improvement for all the members?
          groups_to_check = potential_groups[improved_scores]
          lambdas_to_check = potential_lambdas[improved_scores]
          ds_to_check = potential_groups_ds[improved_scores]
          ids_to_check = groups_ids_to_join[improved_scores]
          better_groups = falses(nb_better_grp)
          #println(nb_better_grp, " ", length(groups_to_check))
          for j = 1:nb_better_grp
            #println("... ", j, "/", nb_better_grp)
            better_groups[j] = check_if_valid_group(p, groups_to_check[j], ds_to_check[j])
          end

          if sum(better_groups)>0
            valid_groups = groups_to_check[better_groups]
            valid_lambdas = lambdas_to_check[better_groups]
            valid_ds = ds_to_check[better_groups]
            valid_ids = ids_to_check[better_groups]
            #join maximum scoring group

            new_group_idx = find(x->x==maximum(valid_ds), valid_ds)[end]#[end]
            #new_group_idx = find(x->x==rand(valid_ds), valid_ds)[1]#[end]

            new_group = valid_groups[new_group_idx]
            new_lambda = valid_lambdas[new_group_idx]
            ids_to_merge_with = valid_ids[new_group_idx]
            agent_ids = [a.id for a in new_group]
            #println("agent ", agent.id, " from group ", target_group_id,
            #        " would benefit of making a group ", agent_ids,
            #        " with agents from groups ", ids_to_merge_with)
            #remove agent from original group ...
            group_remove_agents!(p, target_group_id, agent)
            #println("... agent ", agent.id, " removed from group ", target_group_id)
            #add agents to the biggest groups of the ids to merge with
            grps = filter(g->g.id in ids_to_merge_with, p.groups_list)
            gsizes = [length(g.agents_id_list) for g in grps]
            final_grp_id = grps[find(x->x==maximum(gsizes), gsizes)[1]].id
            group_add_agents!(p, final_grp_id, [agent])
            #println("... agent ", agent.id, " added to group ", final_grp_id)

            #update laziness value of the final group
            final_g = population_get_group(p, final_grp_id)
            final_g.laziness_value = new_lambda

            for gp in ids_to_merge_with
              if gp==final_grp_id
                continue
              else
                g_start = population_get_group(p, gp)
                group_add_agents!(p, final_grp_id, g_start.agents_id_list)
                #println("... agents ", g_start.agents_id_list, " added to group ", final_grp_id)
                g_start.agents_id_list = Set(Int64[]) #clear agent_id_list
                #println("... agents list of group ", g_start.id, " cleared")
              end
            end
            movement+=1
          end
        end

        #if movement occured, break the loop and go back to start
        if movement!=0
          break
        end

      end

    end
    pattern_end = population_get_all_groups_dom_scores(p)
  end
end




function population_random_coalitions_freerider_lossnerve!(p::Population)

  pattern_end = [0. for i in 1:2]
  pattern_start = [1. for i in 1:2]
  while !reached_equilibrium(pattern_end, pattern_start) #loop until agents are no more changing groups
    pattern_start = population_get_all_groups_dom_scores(p)
    ds = copy(pattern_start) #population_get_all_groups_dom_scores(p)
    #non_empty_indices = find(x->x!=0., ds)
    glist = p.groups_list#[non_empty_indices]
    #ds = ds[non_empty_indices]
    current_groups = [(Float64(score), Int64(grp.id), grp) for (score, grp) in zip(ds, glist)] #score/id/group
    sort!(current_groups, by=x->x[1]) #sort by dom score

    nb_groups = length(current_groups)
    groups_agents_and_ids = [(shuffle(filter(a->a.id in g.agents_id_list, p.agents_list)), id) for (s,id,g) in current_groups]

    movement = 0

    #iterate over the sorted list of groups
    for i = 1:nb_groups
      #tic()
      (target_group,target_group_id) = groups_agents_and_ids[i]
      if length(target_group)==0 #if empty group, iterate over next group
        continue
      end
      others_groups_and_ids = deleteat!(copy(groups_agents_and_ids), i)

      #############
      ##get rid of all but 1 empty group
      ##prevent doing several time the same calculation while still
      ##ensuring that agent will leave group if the individual score is higher
      ##than ds of any group he could be part of (since solo group will be checked)
      ############
      unique_combination = []
      pos_to_keep = []
      for m in 1:nb_groups-1
        g2check = others_groups_and_ids[m]
        if g2check[1] in unique_combination
          continue
        else
          push!(pos_to_keep, m)
        end
        push!(unique_combination, g2check[1])
      end
      others_groups_and_ids = [others_groups_and_ids[y] for y in pos_to_keep]

      ###### Note:
      ###### Evaluate option of doing the shuffle of groups to make alliances with
      ###### at each new agent, instead of each new groups

      shuffle!(others_groups_and_ids) #shuffle order of other groups to make alliances with

      #add target group at the end so the full population will also be checked
      #in case one of the agents of the target group is in none of the others groups
      push!(others_groups_and_ids, (target_group,target_group_id))

      others_groups = [g for (g,id) in others_groups_and_ids]
      others_ids = [id for (g,id) in others_groups_and_ids]
      #push!(others_groups, target_group)
      #push!(others_ids, target_group_id)
      others_groups_agents = vcat(others_groups...) #agents from random group order in 1 array

      gsizes = [length(grp) for grp in others_groups]
      sum_gsizes = cumsum(gsizes)

      #list of unique agents to form coalition with and IDs of groups of those agents
      groups_to_join = [Set(others_groups_agents[1:k]) for k in sum_gsizes]
      groups_ids_to_join = [Vector{Int}(others_ids[1:k]) for k in 1:length(gsizes)]
      #println(groups_ids_to_join)
      #each agent of target group tries to defect to others groups or groups of groups
      for agent in target_group
        #current_dom_score of agent
        current_dom_scores = agent_get_current_dom_scores(p, agent)
        current_max_dom_score = maximum(current_dom_scores)

        #scores of the other groups with the agent added
        potential_groups = [push!(copy(grp), agent) for grp in groups_to_join]
        ###### LAZINESS FACTOR for potential groups
        potential_lambdas = rand(p._distribution_groups_laziness, length(potential_groups))
        for l in 1:length(potential_groups)
          if length(potential_groups[l])==1
            potential_lambdas[l]=1.
          end
        end
        #corresponding DS for potential groups
        potential_groups_ds = [get_dom_score_agents_list(list, p, lambda) for (list, lambda) in zip(potential_groups, potential_lambdas)]


        #check if better than current max score
        improved_scores = compare_float_vector(current_max_dom_score, potential_groups_ds) #this is a boolean array
        #improved_groups = potential_groups[improved_scores]

        nb_better_grp = sum(improved_scores)

        if nb_better_grp>0
          ##which of those groups provide improvement for all the members?
          groups_to_check = potential_groups[improved_scores]
          lambdas_to_check = potential_lambdas[improved_scores]
          ds_to_check = potential_groups_ds[improved_scores]
          ids_to_check = groups_ids_to_join[improved_scores]
          better_groups = falses(nb_better_grp)
          #println(nb_better_grp, " ", length(groups_to_check))
          for j = 1:nb_better_grp
            #println("... ", j, "/", nb_better_grp)
            better_groups[j] = check_if_valid_group(p, groups_to_check[j], ds_to_check[j])
          end

          if sum(better_groups)>0
            #50% of the time, loss of nerve
            loss_nerve = rand()
            if loss_nerve<0.5
              break
            end

            valid_groups = groups_to_check[better_groups]
            valid_lambdas = lambdas_to_check[better_groups]
            valid_ds = ds_to_check[better_groups]
            valid_ids = ids_to_check[better_groups]
            #join maximum scoring group

            new_group_idx = find(x->x==maximum(valid_ds), valid_ds)[end]#[end]
            #new_group_idx = find(x->x==rand(valid_ds), valid_ds)[1]#[end]

            new_group = valid_groups[new_group_idx]
            new_lambda = valid_lambdas[new_group_idx]
            ids_to_merge_with = valid_ids[new_group_idx]
            agent_ids = [a.id for a in new_group]
            #println("agent ", agent.id, " from group ", target_group_id,
            #        " would benefit of making a group ", agent_ids,
            #        " with agents from groups ", ids_to_merge_with)
            #remove agent from original group ...
            group_remove_agents!(p, target_group_id, agent)
            #println("... agent ", agent.id, " removed from group ", target_group_id)
            #add agents to the biggest groups of the ids to merge with
            grps = filter(g->g.id in ids_to_merge_with, p.groups_list)
            gsizes = [length(g.agents_id_list) for g in grps]
            final_grp_id = grps[find(x->x==maximum(gsizes), gsizes)[1]].id
            group_add_agents!(p, final_grp_id, [agent])
            #println("... agent ", agent.id, " added to group ", final_grp_id)

            #update laziness value of the final group
            final_g = population_get_group(p, final_grp_id)
            final_g.laziness_value = new_lambda

            for gp in ids_to_merge_with
              if gp==final_grp_id
                continue
              else
                g_start = population_get_group(p, gp)
                group_add_agents!(p, final_grp_id, g_start.agents_id_list)
                #println("... agents ", g_start.agents_id_list, " added to group ", final_grp_id)
                g_start.agents_id_list = Set(Int64[]) #clear agent_id_list
                #println("... agents list of group ", g_start.id, " cleared")
              end
            end
            movement+=1
          end
        end

        #if movement occured, break the loop and go back to start
        if movement!=0
          break
        end

      end

    end
    pattern_end = population_get_all_groups_dom_scores(p)
  end
end

#=
##################
##################
##################
##################
##################
##################
##################
##################

function sample_k!(v::Vector, k::Integer)
  n = length(v)
  for i = 1:k
    j = rand(i:n)
    v[i], v[j] = v[j], v[i]
  end
end

function shuffle_between!(v::Vector, k_start::Integer, k_end::Integer)
  for i = k_start:k_end
    j = rand(i:k_end)
    v[i], v[j] = v[j], v[i]
  end
end

function get_ranked_shuffled_agents_list(p::Population)
  #return (Grp DomScore, Grp ID, shuffled Grp Agent List) of non empty groups in population
  groups_count = length(p.groups_list)
  grpDS_grpID_grpAgents = Array(Tuple{Float64, Int64, Vector{Agent}}, 0)

  for i in 1:groups_count
    g = p.groups_list[i]
    agents_list = shuffle(filter(a->a.id in g.agents_id_list, p.agents_list))
    if length(agents_list)!=0
      g_ds = get_dom_score_agents_list(agents_list, p.g_max, p.g_kd, p.g_h, g.laziness_value)
      push!(grpDS_grpID_grpAgents, (g_ds, g.id, agents_list))
    end
  end
  sort!(grpDS_grpID_grpAgents, by=x->x[1])
  return grpDS_grpID_grpAgents
end

function get_coalition_list(grpDS_grpID_grpAgents::Array{Tuple{Float64, Int64, Vector{Agent}}})
  #list of unique agents to form coalition with and IDs of groups of those agents
  others_scores, others_ids, others_groups = [zip(grpDS_grpID_grpAgents...)...] #unzip
  others_groups_agents = vcat(others_groups...) #agents from random group order in 1 array

  gsizes = [length(grp) for grp in others_groups]
  sum_gsizes = cumsum(gsizes)

  #list of unique agents to form coalition with and IDs of groups of those agents
  groups_to_join = Array(Set{Agent}, length(gsizes))
  groups_ids_to_join = Array(Vector{Int}, length(gsizes))
  for k in 1:length(gsizes)
    groups_to_join[k] = Set(others_groups_agents[1:sum_gsizes[k]])
    groups_ids_to_join[k] = [others_ids[1:k]...]
  end
  return (groups_to_join, groups_ids_to_join)
end

function population_random_coalitions2_freeriders!(p::Population)

  pattern_end = [0. for i in 1:2] #make start and end different so algorithm will take place
  pattern_start = [1. for i in 1:2]

  while !reached_equilibrium(pattern_end, pattern_start) #loop until agents are no more changing groups
    pattern_start = population_get_all_groups_dom_scores(p) #initial state

    grpDS_grpID_grpAgents = get_ranked_shuffled_agents_list(p) #(score/id/agents list)
    nb_groups = length(grpDS_grpID_grpAgents)

    #iterate over the sorted list of groups
    for i = 1:nb_groups
      (target_score, target_group_id, target_group) = grpDS_grpID_grpAgents[i]
      #target_group = grpDS_grpID_grpAgents[i]

      others_groups_and_ids = grpDS_grpID_grpAgents[1:end]#copy(grpDS_grpID_grpAgents)
      #shuffle order of other groups to make alliances with while moving target group
      # at the end so the full population will also be checked (in case one of the agents of the target group is in none of the others group)
      others_groups_and_ids[end], others_groups_and_ids[i] = others_groups_and_ids[i], others_groups_and_ids[end]
      shuffle_between!(others_groups_and_ids, 1, nb_groups-1)

      (groups_to_join, groups_ids_to_join) = get_coalition_list(others_groups_and_ids)

      #each agent of target group tries to defect to others groups or groups of groups
      for agent in target_group#[3]#position 3 is the agents list
        #current_dom_score of agent
        current_dom_scores = agent_get_current_dom_scores(p, agent)
        current_max_dom_score = maximum(current_dom_scores)

        #potential groups to check
        potential_groups = [push!(copy(grp), agent) for grp in groups_to_join]
        ###### LAZINESS FACTOR for potential groups
        potential_lambdas = rand(p._distribution_groups_laziness, length(potential_groups))
        for l in 1:length(potential_groups)
          if length(potential_groups[l])==1
            potential_lambdas[l]=1.
          end
        end
        #corresponding DS for potential groups
        potential_groups_ds = [get_dom_score_agents_list(list, p, lambda) for (list, lambda) in zip(potential_groups, potential_lambdas)]

        #check if better than current max score
        improved_scores = compare_float_vector(current_max_dom_score, potential_groups_ds) #this is a boolean array
        #improved_groups = potential_groups[improved_scores]

        nb_better_grp = sum(improved_scores)

        if nb_better_grp>0 #if valid potential groups exist
          ##which of those groups provide improvement for all the members?
          groups_to_check = potential_groups[improved_scores]
          lambdas_to_check = potential_lambdas[improved_scores]
          ds_to_check = potential_groups_ds[improved_scores]
          ids_to_check = groups_ids_to_join[improved_scores]
          better_groups = falses(nb_better_grp)
          #println(nb_better_grp, " ", length(groups_to_check))
          for j = 1:nb_better_grp
            #println("... ", j, "/", nb_better_grp)
            better_groups[j] = check_if_valid_group(p, groups_to_check[j], ds_to_check[j])
          end

          if sum(better_groups)>0
            valid_groups = groups_to_check[better_groups]
            valid_lambdas = lambdas_to_check[better_groups]
            valid_ds = ds_to_check[better_groups]
            valid_ids = ids_to_check[better_groups]
            #join maximum scoring group

            new_group_idx = find(x->x==maximum(valid_ds), valid_ds)[end]#[end]
            #new_group_idx = find(x->x==rand(valid_ds), valid_ds)[1]#[end]

            new_group = valid_groups[new_group_idx]
            new_lambda = valid_lambdas[new_group_idx]
            ids_to_merge_with = valid_ids[new_group_idx]
            agent_ids = [a.id for a in new_group]
            #println("agent ", agent.id, " from group ", target_group_id,
            #        " would benefit of making a group ", agent_ids,
            #        " with agents from groups ", ids_to_merge_with)
            #remove agent from original group ...
            group_remove_agents!(p, target_group_id, agent)
            #println("... agent ", agent.id, " removed from group ", target_group_id)
            #add agents to the biggest groups of the ids to merge with
            grps = filter(g->g.id in ids_to_merge_with, p.groups_list)
            gsizes = [length(g.agents_id_list) for g in grps]
            final_grp_id = grps[find(x->x==maximum(gsizes), gsizes)[1]].id
            group_add_agents!(p, final_grp_id, [agent])
            #println("... agent ", agent.id, " added to group ", final_grp_id)

            #update laziness value of the group
            final_g = population_get_group(p, final_grp_id)
            final_g.laziness_value = new_lambda

            for gp in ids_to_merge_with
              if gp==final_grp_id
                continue
              else
                g_start = population_get_group(p, gp)
                group_add_agents!(p, final_grp_id, g_start.agents_id_list)
                #println("... agents ", g_start.agents_id_list, " added to group ", final_grp_id)
                g_start.agents_id_list = Set(Int64[]) #clear agent_id_list
                #println("... agents list of group ", g_start.id, " cleared")
              end
            end
          end
        end
      end
      #print(i, ' ')
      #toc()
    end
    pattern_end = population_get_all_groups_dom_scores(p)
  end
end
=#
