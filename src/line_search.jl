function line_search{T<:AbstractFloat}(f::Function,x::Array{T,1},d::Array{T,1})

  alpha = 1.0
  beta  = 0.0001
  tau   = 0.5

  decent = derivative(f,x)*d
  if decent[1] > 0.0
    error("d is not a decent direction")
  end

  while f(x+alpha*d) > f(x) + alpha*beta*decent[1]
    alpha = tau*alpha
  end

  return alpha

end
