function newton_raphson(f::Function,x::Array{Float64,1},tol::Float64,maxiters::Int64)

	x_new = similar(x)

  iter = 0
  len = Inf
  while len > tol

    g = derivative(f,x)
    h = hessian(f,x)

    x_new = x - vec(h\g')

    len = maximum(abs,x_new-x)
    x = copy(x_new)

    iter += 1
    if iter >= maxiters
      break
    end

  end

  return x_new, f(x_new), iter

end
