function pattern_search{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},tol::T,maxiters::S)

  n = length(x)

  const lambda = 0.5

  retcode = true

  # Create initial step-lengths

  dh = zeros(T,n,n)
  for i = 1:n
    dh[i,i] = d[i]
  end

  x_old = x

  iters = 0
  while true

    for i = 1:n

      if f(x+dh[:,i]) < f(x) && f(x+dh[:,i]) < f(x-dh[:,i])
        step = 1
        while f(x+(step+1)*dh[:,i]) < f(x+step*dh[:,i])
          step += 1
        end
        x += step*dh[:,i]
      elseif f(x-dh[:,i]) < f(x) && f(x-dh[:,i]) <= f(x+dh[:,i])
        step = -1
        while f(x+(step-1)*dh[:,i]) < f(x+step*dh[:,i])
          step += -1
        end
        x += step*dh[:,i]
      end

    end

    len = maxabs(x-x_old)

    x_old = x
    dh = lambda*dh

    iters += 1

    if maxabs([len;diag(dh)]) < tol
      break
    end

    if iters == maxiters
      retcode = false
      break
    end

  end

  return x, f(x), retcode

end
