abstract oneDoperator
include("partialDifferentiate/partialD.jl")


abstract chooseDimensionOperator


immutable byDimensionN{TypeOperator <: oneDoperator} <: chooseDimensionOperator
  n :: Int64
  functor :: TypeOperator
end

function resize_dims(n :: Int64, dims::Tuple{Vararg{Int64}})
  (reduce(*, dims[1:pf.n -1]), dims[pf.n], reduce(*, dims[pf.n+1:end])) 
end


function eval(pf :: byDimensionN, data :: Array{Float64})
  dims_data = size(pf)
  dataResize = reshape(data, dims_data)
  res_anexample = eval(pf.functor, collect(dataResize[1,:,1]))
  res = fill(0.0 , size(dataResize,1), length(res_anexample), size(dataResize,3))
  for i1 = 1:size(res,1)
    for i3 = 1:size(res,3)
      res[i1,:,i3] = eval(pf.functor, collect(dataResize[i1,:,i3]))
    end
  end
  reshape(res, (dims_data[1:n-1]..., size(res,2), dims_data[n+1:end]...))
end


