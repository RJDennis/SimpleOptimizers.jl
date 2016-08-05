function select_pair(p)

  i = rand(1:p,2)
  if i[1]==i[2]
    i = select_pair(p)
  end

  return i

end

function shuffle_crossover(n,parent1,parent2)

  parent1 = copy(parent1)
  parent2 = copy(parent2)

  bd = round(Int,rand(n))
  for i = 1:n
    if bd[i] == 1
      parent1[i], parent2[i] = parent2[i], parent1[i]
    end
  end

  return parent1,parent2

end

function arithmetic_crossover(parent1,parent2)

  c = rand()
  p1 = c*parent1     + (1-c)*parent2
  p2 = (1-c)*parent1 + c*parent2

  return p1,p2

end

function single_point_crossover(n,parent1,parent2)

  parent1 = copy(parent1)
  parent2 = copy(parent2)

  cut = rand(1:(n-1))
  parent1[cut+1:n], parent2[cut+1:n] = parent2[cut+1:n], parent1[cut+1:n]

  return parent1,parent2

end

function crossover(phi,n,parent1,parent2)

  cross = rand()
  if cross < phi
    method = rand(1:3)
    if method == 1
      (c1,c2) = shuffle_crossover(n,parent1,parent2)
    elseif method == 2
      (c1,c2) = arithmetic_crossover(parent1,parent2)
    else
      (c1,c2) = single_point_crossover(n,parent1,parent2)
    end
  else
    c1 = parent1
    c2 = parent2
  end

  return c1,c2

end

function elitism(g0,g1,n)

  g1 = copy(g1)

  best_old = indmin(g0[n+1,:])
  best_new = indmin(g1[n+1,:])
  worst_new = indmax(g1[n+1,:])
  if g0[n+1,best_old] < g1[n+1,best_new]
    g1[:,worst_new] = g0[:,best_old]
  end

  return g1

end

function mutation(coef,iters,maxiters)

  mu = 0.15 + 0.33/iters
  if rand() < mu
    b = 2.0
    r1 = rand()
    r2 = rand()
    s  = randn()
    if r1 > 0.5
      coef = coef + s*(1.0-(r2^((1.0-iters/maxiters)^b)))
    else
      coef = coef - s*(1.0-(r2^((1.0-iters/maxiters)^b)))
    end
  end

  return coef

end

function genetic_search(f,x,d,m,tol,maxiters)

  x = copy(x)
  n = length(x)

  phi = 0.95
  iters = 0
  retcode = true

  # Ensure the population size is even and non-empty

  if isodd(m)
    m += 1
  end

  # Construct the initial parent population

  domain_lower = x - d
  domain_upper = x + d

  g0 = domain_lower .+ rand(n,m).*(domain_upper - domain_lower)
  g0 = vcat(g0,zeros(1,m))

  # Evaluate and store the performance of each parent

  for i = 1:m
    g0[n+1,i] = f(g0[1:n,i])
  end

  # Initialize the child population

  while true

    g1 = zeros(n+1,m)

    for i = 1:round(Int,m/2)

      pair1 = g0[:,select_pair(m)]
      pair2 = g0[:,select_pair(m)]

      parent1 = pair1[:,indmin(pair1[n+1,:])]
      parent2 = pair2[:,indmin(pair2[n+1,:])]

      (child1,child2) = crossover(phi,n,parent1,parent2)

      # Mutation

      for j = 1:n
        child1[j] = mutation(child1[j],iters,maxiters)
        child2[j] = mutation(child2[j],iters,maxiters)
      end

      child1[n+1] = f(child1[1:n])
      child2[n+1] = f(child2[1:n])

      family = [parent1 parent2 child1 child2]
      sort_index = sortperm(collect(family[n+1,:]))
      sorted_family = family[:,sort_index]
      g1[:,2*(i-1)+1:2*i] = sorted_family[:,3:4]

    end

    g1 = elitism(g0,g1,n)

    iters += 1
    if maxabs(mean(g0,2)-mean(g1,2)) < tol
      g0 = copy(g1)
      break
    end

    g0 = copy(g1)

    if iters == maxiters
      retcode = false
      break
    end

  end

  x = g0[1:n,indmin(g0[n+1,:])]
  f = g0[n+1,indmin(g0[n+1,:])]

  return x,f,retcode

end
