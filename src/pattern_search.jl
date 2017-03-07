function pattern_search{T<:AbstractFloat}(f::Function,x::Array{T,1},h::Array{T,1},tol::T,maxiters::Integer)

  x_new = similar(x)
  f0    = f(x)

  lambda = 0.5

  iters = 0
  len = Inf
  while len > tol

    for i = 1:length(x)

      x_new[i] = x[i] + h[i]
      if f(x_new) < f0
        f0 = f(x_new)
      else
        x_new[i] = x[i] - h[i]
        if f(x_new) < f0
          f0 = f(x_new)
        else
          x_new[i] = x[i]
        end
      end

    end

    if x_new == x
      h = lambda*h
    else
      d = x_new-x
      alpha = line_search_simple(f,x,d)
      x_new = x + alpha*d
    end

    len = maximum(abs,h)
    x = copy(x_new)

    iters += 1
    if iters >= maxiters
      break
    end

  end

  return x, f(x), iters

end
