function initL{Tboundary <:BoundaryStrategy }(
  pf :: PartialD_OP{PartialD{1}, StencilsnN{3,3}, Tboundary}, 
  data :: Array{Float64,1}
  )
  a = zeros(data);
  b = zeros(data);
  c = zeros(data);
  a[2:end-1] = 2/3.0
  c[2:end-1] = 1/6.0
  b[2:end-1] = 1/6.0
  a[1] = 1.0; a[end] = 1.0;
  c[1] = 3.0; b[end] = 3.0;
  return(a,b,c)
end

function initL{Tboundary <:BoundaryStrategy }(
  pf :: PartialD_OP{PartialD{2}, StencilsnN{3,3}, Tboundary}, 
  data :: Array{Float64,1}
  )
  a = zeros(data);
  b = zeros(data);
  c = zeros(data);
  a[2:end-1] = 5/6.0
  c[2:end-1] = 1/12.0
  b[2:end-1] = 1/12.0
  a[1] = 1.0; a[end] = 1.0;
  c[1] = 11.0; b[end] = 11.0;
  return(a,b,c)
end

function initR{Tboundary <:BoundaryStrategy }(
  pf :: PartialD_OP{PartialD{1}, StencilsnN{3,3}, Tboundary}, 
  data :: Array{Float64,1}
  )
  res = zeros(data)
  res[2:end-1] = (data[3:end] - data[1:end-2])/2
  pboundary = [-17/6, 3/2 , 3/2, -1/6]
  res[1] = sum( pboundary .* data[1:4])
  res[end] =  - sum( pboundary .* data[end:-1:end-3])
  return res
end

function initR{Tboundary <:BoundaryStrategy }(
  pf :: PartialD_OP{PartialD{2}, StencilsnN{3,3}, Tboundary}, 
  data :: Array{Float64,1}
  )
  res = zeros(data)
  res[2:end-1] = data[3:end] + data[1:end-2] - 2*data[2:end-1]
  pboundary = [13, -27 , 15, -1]
  res[1] = sum( pboundary .* data[1:4])
  res[end] =   sum( pboundary .* data[end:-1:end-3])
  return res
end

