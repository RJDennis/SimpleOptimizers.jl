function construct_initial_simplex{T<:AbstractFloat}(x::Array{T,1},step_size::T)

  n   = length(x)
  p   = Array{T}(n,n+1)
  inc = Array{T}(n)

  for i = 1:n
	  if x[i] == 0.0
	    inc[i] = step_size
		else
		  inc[i] = step_size*x[i]
		end
	end

	p[:,1:(n+1)] = x.*ones(n,n+1)
	p[:,2:(n+1)] += inc.*eye(n)

  return p

end

function compute_reflection_point{T<:AbstractFloat,S<:Integer}(f::Function,y::Array{T,1},p::Array{T,2},n::S,alpha::T)

  index_table = sortperm(y)
	index_worst = index_table[n+1]

  pbar = (sum(p,2) - p[:,index_worst])/n
  pr   = (1.0+alpha)*pbar - alpha*p[:,index_worst]
  yr  = f(pr)

  return index_table, pbar, pr, yr

end

function nelder_mead{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},step_size::T,tol::T,maxiters::S)

  # Some initializations

  alpha = 1.0
  beta  = 0.5
  gam   = 2.0

  n       = length(x)
	x_new   = Array{T}(n)
	f_opt   = 0.0

	# Construct the initial simplex

  p = construct_initial_simplex(x,step_size)

	# Evaluate the function at each vertex

  y = Array{T}(n+1)

  for i = 1:(n+1)
	  y[i] = f(p[:,i])
	end

  # Begin the algorithm

  iters = 0
	retcode = true

  while true

    (index_table,pbar,pr,yr) = compute_reflection_point(f,y,p,n,alpha)

		index_best         = index_table[1]
  	index_second_worst = index_table[n]
	  index_worst        = index_table[n+1]

    if yr < y[index_second_worst] && yr > y[index_best] # Do reflection
      p[:,index_worst] = pr
      y[index_worst]   = yr
    elseif yr < y[index_best]                           # Do expansion
      pe = pbar + gam*(pbar-p[:,index_worst])
      ye = f(pe)
      if ye < y[index_best]
        p[:,index_worst] = pe
				y[index_worst]   = ye
			else
        p[:,index_worst] = pr
        y[index_worst]   = yr
      end
    else                                                # Compute best contraction point
      pc1 = pbar-beta*(pbar-p[:,index_worst])
      pc2 = pbar+beta*(pbar-p[:,index_worst])
      if f(pc1) < f(pc2)
        pc = copy(pc1)
      else
        pc = copy(pc2)
      end
      yc = f(pc)
		  if yc < y[index_worst]                            # Do contraction
        p[:,index_worst] = pc
        y[index_worst]   = yc
      else                                              # Do reduction
        for i = 2:(n+1)
          p[:,index_table[i]] = p[:,index_best]+beta*(p[:,index_table[i]]-p[:,index_best])
        end
        for i = 1:(n+1)
	        y[i] = f(p[:,i])
	      end
      end
    end

    len = 0.0
    mean_p = mean(p,2)
	  for i = 1:n
	    for j = 1:(n+1)
		    len = maximum(abs,[len,(p[i,j]-mean_p[i])])
		  end
	  end

    if len < tol
	    x_new  = vec(mean_p)  # vec() is to convert x to 1d array
	    f_opt = mean(y)
		  break
    end

    iters += 1
    if iters >= maxiters
	    x_new   = vec(mean_p)  # vec() is to convert x to 1d array
	    f_opt   = mean(y)
		  break
		end

  end

	return x_new, f_opt, iters

end
