function hjexplore{T<:AbstractFloat}(f::Function,x0::Array{T,1},d::Array{T,1})

  n = length(x0)

  x = copy(x0)
  f0 = f(x0)

  for i = 1:n

    x[i] += d[i]
    ff = f(x)

    if ff < f0

      f0 = ff

    else

      x[i] = x0[i]-d[i]
      ff = f(x)

      if ff < f0

        f0 = ff

      else

        x[i] = x0[i]

      end

    end

  end

  return x

end

function hooke_jeeves{T<:AbstractFloat,S<:Integer}(f::Function,x0::Array{T,1},d::Array{T,1},tol::T,maxiters::S)

  x = copy(x0)

  retcode = true
  iters = 0

  while maxabs(d) > tol

    y = hjexplore(f,x,d)

    if maxabs(y-x) < tol

      d = 0.5*d

    else

      x = 2*y-x
      z = hjexplore(f,x,d)

      if maxabs(z-x) < tol

        x = y

      else

        x = z

      end

    end

    iters += 1

    if iters == maxiters

      retcode = false
      break

    end

  end

  return x, f(x), retcode

end
