function hjexplore{T<:AbstractFloat}(f::Function,x::Array{T,1},d::Array{T,1})

  n = length(x)
  x_new = copy(x)

  f0 = f(x)
  for i = 1:n
    x_new[i] += d[i]
    ff = f(x_new)
    if ff < f0
      f0 = ff
    else
      x_new[i] = x[i]-d[i]
      ff = f(x_new)
      if ff < f0
        f0 = ff
      else
        x_new[i] = x[i]
      end
    end
  end

  return x_new

end

function hooke_jeeves{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},tol::T,maxiters::S)

  x_new = copy(x)

  iters = 0
  len = Inf
  while len > tol

    y = hjexplore(f,x_new,d)
    if maximum(abs,y-x_new) < tol
      d = 0.5*d
    else
      x_new = 2*y-x_new
      z = hjexplore(f,x_new,d)
      if maximum(abs,z-x_new) < tol
        x_new = copy(y)
      else
        x_new = copy(z)
      end
    end

    len = maximum(abs,d)

    iters += 1
    if iters >= maxiters
      break
    end

  end

  return x_new, f(x_new), iters

end
