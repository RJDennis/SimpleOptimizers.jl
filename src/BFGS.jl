function bfgs{T<:AbstractFloat}(f::Function,x::Array{T,1},tol::T,maxiters::Integer)

  n = length(x)
	x_new = similar(x)

  iters = 0
  len = Inf

  g = derivative(f,x)
  h = inv(hessian(f,x))

  while len > tol

	  xi = -vec(h*g')

		for i = 1:n
		  x_new[i] = x[i] + xi[i]
		end

    gn = derivative(f,x_new)
    dg = (gn - g)'
    hdg = h*dg
    fac = dg'xi
    fae = dg'hdg
    if (fac != 0.0) && (fae != 0.0)
      hn = h + (1./fac).*xi*xi' - (1./fae).*hdg*hdg' + fae.*((1./fac).*xi - (1./fae).*hdg).*((1./fac).*xi - (1./fae).*hdg)'
    else
      break
    end

    len = maximum(abs,x_new-x)
    x = copy(x_new)
    g = copy(gn)
    h = copy(hn)

    iters += 1
    if iters >= maxiters
      break
    end

  end

  return x, f(x), iters

end
