function conjugate_gradient{T<:AbstractFloat}(f::Function,b0::Array{T,1},tol::T,maxiters::Integer,method::AbstractString="FR")

  iters   = 0
  retcode = true

  const beta   = 0.0001
	const sigma  = 0.40
	const rho    = 0.9
	const alpha0 = 1.0

  b0_n = similar(b0)
  fx1 = Float64

	len = Inf
	while true

 	   g = collect(derivative(f,b0)') # collect() is used to convert to 1d array
	   s = -g
		 fx = f(b0)

		# perform line search

		alpha = alpha0
		count = 1
		while count <= 20

			b0_n = b0 + alpha*s
			fx1 = f(b0_n)
			line_count = 1
		  if fx1 > fx + alpha*beta*sum(g's) # check first Wolfe condition # sum is used to convert single-element array to float
				line_count += 1
				if line_count < 20
					alpha = alpha*rho
				else
					alpha = -alpha0
				end
			else
				g1 = collect(derivative(f,b0_n)') # collect() is used to convert to 1d array
			  if maxabs(g1's) > sigma*maxabs(g's) # check second Wolfe condition # maxabs is used instead of abs to convert single-element arrays to floats
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
            elseif mrthod == "HS"
              delta = maximum([(g1'g1-g1'g)/(s'g1-s'g);0.0]) # Hestenes-Stiefel
            elseif method == "DY"
              delta = maximum([(g1'g1)/(s'g1-s'g);0.0])      # Dai–Yuan
            else
              println("conjugate_gradient: Invalid method")
              throw()
            end
					end
          s = -g1 + delta*s
					len = maxabs(g1)
					if len < tol
						break
					end
					g = g1
				  b0 = b0_n
					count += 1
				end
			end

		end

    iters += 1

	  if len < tol
	    break
	  end

    if iters == maxiters
		  retcode = false
      break
	  end

	end

  return b0_n, fx1, retcode

end
