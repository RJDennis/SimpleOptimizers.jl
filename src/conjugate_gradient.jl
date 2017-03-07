function conjugate_gradient{T<:AbstractFloat}(f::Function,x::Array{T,1},tol::T,maxiters::Integer;method::AbstractString="FR")

  x_new = similar(x)

  iters   = 0
	len = Inf
	while len > tol

 	  g = vec(derivative(f,x)') # vec() is used to convert to 1d array
	  s = -g

    alpha = line_search_wolfe(f,x,s)
    x_new = x + alpha*s
    g1 = vec(derivative(f,x_new)') # vec() is used to convert to 1d array

		if sum(g1'g/g1'g1) > 0.1  # perform restart # sum is used to convert single-element array to float
			delta = 0.0
		else
      if method == "FR"
        delta = maximum(g1'g1/g'g)                     # Fletcher-Reeves
      elseif method == "PR"
        delta = maximum([(g1'g1-g1'g)/g'g;0.0])        # Polak-Ribiere
      elseif method == "HS"
        delta = maximum([(g1'g1-g1'g)/(s'g1-s'g);0.0]) # Hestenes-Stiefel
      elseif method == "DY"
        delta = maximum([(g1'g1)/(s'g1-s'g);0.0])      # Dai–Yuan
      else
        throw("conjugate_gradient: Invalid method")
      end
		end

    s = -g1 + delta*s

    len = maximum(abs,g1)
    g = copy(g1)
	  x = copy(x_new)

    iters += 1
    if iters >= maxiters
      break
    end

	end

  return x, f(x), iters

end
