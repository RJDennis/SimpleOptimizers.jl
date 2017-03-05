function conjugate_gradient{T<:AbstractFloat}(f::Function,x::Array{T,1},tol::T,maxiters::Integer;method::AbstractString="FR")

  beta   = 0.0001
	sigma  = 0.40
	rho    = 0.9
	alpha0 = 1.0

  x_new = similar(x)
  fx1 = Float64

  iters   = 0
	len = Inf
	while true

 	  g = vec(derivative(f,x)') # vec() is used to convert to 1d array
	  s = -g
		fx = f(x)

		# perform line search

		alpha = alpha0
		count = 1
    len_inner = Inf
		while count <= 20 && len_inner > tol

			x_new = x + alpha*s
			fx1 = f(x_new)
			line_count = 1
		  if fx1 > fx + alpha*beta*sum(g's) # check first Wolfe condition # sum is used to convert single-element array to float
				line_count += 1
				if line_count < 20
					alpha = alpha*rho
				else
					alpha = -alpha0
				end
			else
				g1 = vec(derivative(f,x_new)') # vec() is used to convert to 1d array
			  if maximum(abs,g1's) > sigma*maximum(abs,g's) # check second Wolfe condition # maxabs is used instead of abs to convert single-element arrays to floats
					line_count += 1
					if line_count < 20
					  alpha = alpha*rho
				  else
					  alpha = -alpha0
				  end
				else
					# line search complete
				  # update optimization point
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
              println("conjugate_gradient: Invalid method")
              throw()
            end
					end
          s = -g1 + delta*s
					len = maximum(abs,g1)
          if len < tol
            break
          end
					g = copy(g1)
				  x = copy(x_new)
					count += 1
				end
			end

		end

	  if len < tol
	    break
	  end

    iters += 1
    if iters >= maxiters
      break
	  end

	end

  return x, fx1, iters

end
