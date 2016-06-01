type SOR_onestep{T <: Localized_Functor} <: GridMatrixIter_one
  A :: T
  OL :: Tuple{Vararg{Int64}}
end

function eval!{T<:Localized_Functor,n}(
  Ja :: SOR_onestep{T}
  , x :: Array{Float64,n}
  , y :: Array{Float64,n}
  , param :: Tuple{Vararg{Array{Float64,n}}}
  , global_param :: Tuple = ()
  ; alpha = 1.0, iseq :: Array{Int64,1} = Array(Int64,0)
  , delta_x :: Float64 = 0.5)

  OL = Ja.OL
  irank = map((1:n...)) do ii
    size(x,ii) - 2OL[ii]
  end

  Ais =  get_all_Ai(
           Ja.OL, Ja.A, x, param, global_param; 
           iseq = iseq, delta_x = delta_x, OL_bound = Ja.OL
           )


  updatepoint_range = ( 0 == length(iseq)) ? collect(1:reduce(*,irank)) : iseq
  for (ilogger,ipoint) = enumerate(updatepoint_range)
    iixy = ind2sub(irank,ipoint)

    iixyrange = map((1:n...)) do ii 
      OL[ii]+ ((iixy[ii]-OL[ii]):(iixy[ii]+OL[ii]))
    end

    param_local = map(param) do xx
      xx[iixyrange...]
    end

    xbackground = x[iixyrange...]
    idx = map(x->x+1, OL)

    Ax = eval(Ja.A, xbackground, param_local, idx, global_param)
    Ai = Ais[ilogger]

    y_localerror = y[map(ii -> OL[ii]+iixy[ii], (1:n...))... ] - Ax
    xbackground_idx = map(ii -> OL[ii] + iixy[ii], (1:n...))
    x[xbackground_idx...] += alpha*y_localerror / (Ai+ 0)
  end
  return x
end

include("SOR.jl")
