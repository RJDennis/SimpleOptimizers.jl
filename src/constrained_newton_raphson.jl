function constrained_newton_raphson{T<:AbstractFloat}(f::Function,b0::Array{T,1},b0_lower::Array{T,1},b0_upper::Array{T,1},tol::T,maxiters::Integer)

  iters   = 0
  retcode = true
  n = length(b0)
  b0_n = similar(b0)

  while true

    g = derivative(f,b0)
    h = hessian(f,b0)
    update = h\g'

    for i = 1:n
	    b0_n[i] = b0[i] - update[i]
      if b0_n[i] > b0_upper[i]
        b0_n[i] = b0_upper[i]
      end
      if b0_n[i] < b0_lower[i]
        b0_n[i] = b0_lower[i]
      end
    end

    if maxabs([(b0-b0_n);g']) < tol
      break
    end

    if iters == maxiters
      retcode = false
      break
    end

    b0 = b0_n
    iters += 1

  end

  return b0, f(b0), retcode

end