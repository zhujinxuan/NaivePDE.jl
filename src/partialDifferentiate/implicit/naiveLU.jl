function naiveLU( a:: Array{Float64,1}, b :: Array{Float64,1}, c :: Array{Float64,1})
  #= @assert (size(a) == size(b)) =#
  #= @assert (size(a) == size(c)) =#
  l = zeros(a)
  d = zeros(a)
  u = zeros(a)
  d[1] = a[1]
  u[1:end-1] = c[1:end-1]
  
  for ii = 2:(length(l))
    l[ii] = b[ii]/d[ii-1]
    d[ii] = a[ii] - l[ii]u[ii-1]
  end
  return (l,d,u)
end
export naiveLU

function ForwardThomas( l :: Array{Float64,1}, dataR :: Array{Float64,1})
  y = zeros(dataR)
  y[1] = dataR[1]
  for ii = 2:length(dataR)
    y[ii] = dataR[ii] - l[ii]*y[ii-1]
  end
  return y
end

function BackwardThomas( d :: Array{Float64, 1}, u :: Array{Float64,1}, dataR :: Array{Float64,1})
  x = zeros(dataR)
  x[end] = dataR[end]/d[end]
  for ii = (length(x)-1):-1:1
    x[ii] = (dataR[ii] - u[ii]*x[ii+1])/d[ii]
  end
  return x
end
