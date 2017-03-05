function mesh_search{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},step::T,m::S,tol::T,maxiters::S)

  x = copy(x)
  n = length(x)
  d = Array{T}(n)

  for i = 1:n
	  if x[i] == 0.0
	    d[i] = step
		else
		  d[i] = step*x[i]
		end
	end

  retcode = true
  iters = 0

  while true

    # Create the mesh

    mesh_lower = x-d
    mesh_upper = x+d

    mesh_edges = zeros(m+1,n)
    for i = 1:n
      mesh_edges[:,i] = collect(linspace(mesh_lower[i],mesh_upper[i],m+1))
    end

    mesh_mid_points = ((mesh_edges[1:m,:]+mesh_edges[2:(m+1),:])/2)
    mesh_size = tuple(fill(m,n)...)

    f_min = Inf
    f_min_index = Tuple[]
    x_new = zeros(n)

    for i = 1:m^n
      index_vector = ind2sub(mesh_size,i)
      point = [mesh_mid_points[index_vector[1],1]]
      for j = 2:n
        point = push!(point,mesh_mid_points[index_vector[j],j])
      end
      f_new = f(point)
      if f_new < f_min
        f_min = f_new
        x_new[:] = point
        f_min_index = index_vector
      end
    end

    # Determine if a minimum has been found

    if maximum(abs,x_new-x) < tol
      break
    end
    x = x_new

    # Determine if the iteration limit has been exceeded

    iters += 1
    if iters == maxiters
      retcode = false
      break
    end

    # Construct the new d vector

    for j = 1:n
      d[j] = mesh_edges[f_min_index[j],j]-x[j]
    end

  end

  return x, f(x), iters

end

function mesh_search{T<:AbstractFloat,S<:Integer}(f::Function,x::Array{T,1},d::Array{T,1},m::S,tol::T,maxiters::S)

  x = copy(x)
  d = copy(d)

  n = length(x)
  retcode = true
  iters = 0

  while true

    # Create the mesh

    mesh_lower = x-d
    mesh_upper = x+d

    mesh_edges = zeros(m+1,n)
    for i = 1:n
      mesh_edges[:,i] = collect(linspace(mesh_lower[i],mesh_upper[i],m+1))
    end

    mesh_mid_points = ((mesh_edges[1:m,:]+mesh_edges[2:(m+1),:])/2)
    mesh_size = tuple(fill(m,n)...)

    f_min = Inf
    f_min_index = Tuple[]
    x_new = zeros(n)

    for i = 1:m^n
      index_vector = ind2sub(mesh_size,i)
      point = [mesh_mid_points[index_vector[1],1]]
      for j = 2:n
        point = push!(point,mesh_mid_points[index_vector[j],j])
      end
      f_new = f(point)
      if f_new < f_min
        f_min = f_new
        x_new[:] = point
        f_min_index = index_vector
      end
    end

    # Determine if a minimum has been found

    if maximum(abs,x_new-x) < tol
      break
    end
    x = x_new

    # Determine if the iteration limit has been exceeded

    if iters == maxiters
      retcode = false
      break
    end
    iters += 1

    # Construct the new d vector

    for j = 1:n
      d[j] = mesh_edges[f_min_index[j],j]-x[j]
    end

  end

  return x, f(x), iters

end
