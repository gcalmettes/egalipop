using Distributions
include("model-objects.jl")

####################
# Distributions for alpha and beta
####################

function dist_poisson(n)
  lower, upper = -2, 0
  mu, sigma = -1, 0.5
  dist = TruncatedNormal(mu, sigma, lower, upper)
  param = rand(dist, n)
  for i = 1:n
      param[i] = 0.5*10^param[i]
  end
  return param
end

function dist_normal(n)
  mu, sigma = 0.25, 0.06
  dist = TruncatedNormal(mu, sigma, 0, 1)
  param = rand(dist, n)
  return param
end

function dist_uniform(n)
  lower, upper = 0.005, 0.5
  dist = Uniform(lower, upper)
  param = rand(dist, n)
  return param
end


###########
#genral Functions
##########
function get_dom_score_agents_list(v::Vector{Agent}, g_max::Float64, g_kd::Float64, g_h::Float64, lambda::Float64)
  if length(v)==1
    g=1
  else
    g = g_max/(1+length(v)/g_kd)^g_h
  end
  ds=0.
  for a in v
    ds+=a.dom_score
  end
  ds*=g
  ds*=lambda
  return ds
end
function get_dom_score_agents_list(v::Set{Agent}, g_max::Float64, g_kd::Float64, g_h::Float64, lambda::Float64)
  if length(v)==1
    g=1
  else
    g = g_max/(1+length(v)/g_kd)^g_h
  end
  ds=0.
  for a in v
    ds+=a.dom_score
  end
  ds*=g
  ds*=lambda
  return ds
end
function get_dom_score_agents_list(v::Vector{Agent}, p::Population, lambda::Float64)
  if length(v)==1
    g=1
  else
    g = p.g_max/(1+length(v)/p.g_kd)^p.g_h
  end
  ds=0.
  for a in v
    ds+=a.dom_score
  end
  ds*=g
  ds*=lambda
  return ds
end
function get_dom_score_agents_list(v::Set{Agent}, p::Population, lambda::Float64)
  if length(v)==1
    g=1
  else
    g = p.g_max/(1+length(v)/p.g_kd)^p.g_h
  end
  ds=0.
  for a in v
    ds+=a.dom_score
  end
  ds*=g
  ds*=lambda
  return ds
end

function get_agent_ids(v::Vector{Agent})
  idx = [a.id for a in v]
  return idx
end

function get_agents(a::Vector{Agent}, i::Vector{Int})
  agents = filter(x->x.id in i, a)
  return agents
end

function calcul_gamma(size::Int, g_max::Float64, g_kd::Float64, g_h::Float64)
  if size==1
    gamma = 1
  else
    gamma = g_max/(1+size/g_kd)^g_h
  end
  return gamma
end

function compare_vectors(v1::Vector, v2::Vector ; comparison= <)
  length(v1)==length(v2) || error("Cannot compare objects with different lengths")
  n = length(v1)
  res = falses(n)
  for i in 1:n
    res[i] = comparison(v1[i], v2[i])
  end
  return res
end

function compare_float_vector(s::Float64, ps::Vector{Float64})
    "return boolean result of vector.>float"
    size = length(ps)
    r = falses(size)
    for i = 1:size
        r[i] = ps[i]>s
    end
    return r
end

function check_if_valid_group(p::Population, alist::Vector{Agent}, s::Float64)
    n = length(alist)
    n==0 && return false
    agents_ds = [maximum(agent_get_current_dom_scores(p, a)) for a in alist]
    c = 0
    for i = 1:n
        c += s>=agents_ds[i]
    end
    return c==n
end
function check_if_valid_group(p::Population, alist::Set{Agent}, s::Float64)
    n = length(alist)
    n==0 && return false
    agents_ds = [maximum(agent_get_current_dom_scores(p, a)) for a in alist]
    c = 0
    for i = 1:n
        c += s>=agents_ds[i]
    end
    return c==n
end
