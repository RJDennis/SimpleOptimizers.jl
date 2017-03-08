function luus_jaakola{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},seed::S,tol::T,maxiters::S)

  # d is the radius of the initial search area about the initial x
  # maxiters should be large, especially if d is chosen to be too small

  srand(seed)

  lambda = 0.95
  d_store = d
  x_old = copy(x)

  n = length(x)

  len_outer = Inf
  iters = 0
  while len_outer > tol
    len_inner = Inf
    d = d_store
    while len_inner > tol

      z = rand(n)
      y = x - d + 2*z.*d

      if f(y) < f(x)
        len_inner = maximum(abs,[x-y; f(y)-f(x)])
        x = y
      else
        d = lambda*d
      end

      iters += 1
      if iters >= maxiters
        break
      end

    end

    len_outer = maximum(abs,x-x_old)
    x_old = copy(x)

    if iters >= maxiters
      break
    end

  end

  return x, f(x), iters

end
