function constrained_newton_raphson{T<:AbstractFloat}(f::Function,x::Array{T,1},x_lower::Array{T,1},x_upper::Array{T,1},tol::T,maxiters::Integer)

  n = length(b0)
  x_new = similar(x)

  iters = 0
  len   = Inf

  while len > tol

    g = derivative(f,x)
    h = hessian(f,x)
    update = vec(h\g')

    for i = 1:n
	    x_new[i] = x[i] - update[i]
      if x_new[i] > x_upper[i]
        x_new[i] = x_upper[i]
      end
      if x_new[i] < x_lower[i]
        x_new[i] = x_lower[i]
      end
    end

    len = maximum(abs,x_new-x)
    x = copy(x_new)

    iters += 1
    if iters >= maxiters
      break
    end

  end

  return x, f(x), iters

end
