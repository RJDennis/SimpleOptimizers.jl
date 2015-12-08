function constrained_bfgs{T<:AbstractFloat}(f::Function,b0::Array{T,1},b0_lower::Array{T,1},b0_upper::Array{T,1},tol::T,maxiters::Integer)

  iters = 0
  retcode = true
  n = length(b0)
  b0_n = similar(b0)

  g = derivative(f,b0)
  h = inv(hessian(f,b0))

  while true

	  xi = -h*g'

	  for i = 1:n
		  b0_n[i] = b0[i] + xi[i]
      if b0_n[i] > b0_upper[i]
        b0_n[i] = b0_upper[i]
      end
      if b0_n[i] < b0_lower[i]
        b0_n[i] = b0_lower[i]
      end
    end

    gn = derivative(f,b0_n)
    dg = (gn - g)'
    hdg = h*dg
    fac = dg'xi
    fae = dg'hdg
    if (fac != 0.0) && (fae != 0.0)
      hn = h + (1./fac).*xi*xi' - (1./fae).*hdg*hdg' + fae.*((1./fac).*xi - (1./fae).*hdg).*((1./fac).*xi - (1./fae).*hdg)'
    else
      retcode = false
      break
    end

    iters += 1

    if maxabs([(b0_n-b0);g']) < tol
		  break
		end

    if iters == maxiters
      retcode = false
      break
    end

		g = gn
    h = hn
    b0 = b0_n

  end

  return b0, f(b0), retcode

end
