function quadratic_interpolation{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},tol::T,maxiters::S)

  n = length(x)

  const lambda = 0.4

  x0 = x
  x0_n = Array(T,n)

  retcode = true

  iters = 0
  while true

    h = d.*eye(n)

    for i = 1:n

      while true

        x1 = x0 + h[:,i]
        x2 = x0 + 2*h[:,i]

        f0 = f(x0)
        f1 = f(x1)
        f2 = f(x2)

        if f2 < f1 && f1 < f0 # Move further in this direction

         h[i,i] = 2*h[i,i]

        elseif f2 > f1 && f1 > f0 # Move in the opposite direction and decrease step in case at the minimum already

          h[i,i] = -0.5*h[i,i]

        elseif f2 < f1 && f1 > f0 # Have gone too far, or need to jump over a local maximum

          if f2 < f0
            x0_n[i] = x2
          else
            h[i,i] = lambda*h[i,i]
          end

        else # Ideal situation, update x[i]

          x0_n[i] = x0[i]+h[i,i]*(4*f1-3*f0-f2)/(4*f1-2*f0-2*f2)
          break

        end

      end

    end

    iters += 1

    len = maxabs(x0_n-x0)
    x0 = x0_n

    if len < tol
      break
    end

    if iters == maxiters
      retcode = false
      break
    end

  end

  return x0, f(x0), retcode

end
