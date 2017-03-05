function luus_jaakola{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},range::Array{T,2},tol::T,maxiters::S)

  lambda = 0.95

  n = length(x)

  len = Inf
  iters   = 0
  while true

    for i = 1:n
      d[i] = minimum([d[i],abs(x[i]-range[1,i]),abs(x[i]-range[2,i])])
      d[i] = maximum([d[i],tol])
    end

    z = rand(n)
    y = x - d + 2*z.*d

    if f(y) < f(x)
      len = maximum(abs,[x-y; f(y)-f(x)])
      x = y
    else
      d = lambda*d
    end

    if len < tol
      break
    end

    iters += 1
    if iters >= maxiters
      break
    end

  end

  return x, f(x), iters

end
