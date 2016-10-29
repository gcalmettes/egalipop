function random_split(pool::Vector)
  initial_size=length(pool)
  remaining_size=length(pool)

  indices=shuffle(collect(1:initial_size))

  groups = []
  while remaining_size>0
    rdm = rand(1:remaining_size)
    step=initial_size-remaining_size
    push!(groups, pool[indices[step+1:step+rdm]])
    remaining_size -= rdm
  end
  return groups
end
