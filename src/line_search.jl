function line_search_simple{T<:AbstractFloat}(f::Function,x::Array{T,1},d::Array{T,1})

  alpha = 1.0
  beta  = 0.0001
  tau   = 0.9

  decent = derivative(f,x)*d
  if decent[1] > 0.0
    error("d is not a decent direction")
  end

  while f(x+alpha*d) > f(x) + alpha*beta*decent[1] # check first Wolfe condition
    alpha = tau*alpha
  end

  return alpha

end

function line_search_wolfe{T<:AbstractFloat}(f::Function,x::Array{T,1},d::Array{T,1})

  beta   = 0.0001
	sigma  = 0.40
	tau    = 0.9
	alpha = 1.0

  decent = derivative(f,x)*d
  if decent[1] > 0.0
    error("d is not a decent direction")
  end

  # Check Wolfe conditions

  while f(x+alpha*d) > f(x) + alpha*beta*decent[1] && maximum(abs,derivative(f,x+alpha*d)*d) > sigma*maximum(abs,decent)
    alpha = tau*alpha
  end

  return alpha

end
