function get_all_Ai{T<: Localized_Functor,n}(
  OL :: Tuple{Vararg{Int64}}
  ,p :: T
  , x :: Array{Float64,n}
  , param :: Tuple{Vararg{Array{Float64,n}}}
  , global_param :: Tuple = ()
  ; iseq :: Array{Int64,1} = Array(Int64,0)
  , delta_x :: Float64 = 0.5
  , OL_bound :: Tuple{Vararg{Int64}} = OL)

  
  irank = map((1:n...)) do ii
    size(x,ii) - 2OL_bound[ii]
  end

  updatepoint_range = ( 0 == length(iseq)) ? collect(1:reduce(*,irank)) : iseq

  Ai = fill(0.0,length(updatepoint_range))

  for (ilogger, ipoint) = enumerate(updatepoint_range)
    iixy = ind2sub(irank,ipoint)
    iixyrange = map((1:n...)) do ii 
      OL_bound[ii] + ((iixy[ii]-OL[ii]):(iixy[ii]+OL[ii]))
    end

    param_local = map(param) do xx
      xx[iixyrange...]
    end

    xbackground = x[iixyrange...]
    idx = map(x->x+1, OL)
    

    Ai[ilogger] = begin
      x1background = fill(0.0,size(xbackground)); 
      x1background[idx...] += delta_x
      Ax1 = eval(p, x1background, param_local, idx, global_param)
      x0background = fill(0.0,size(xbackground)); 
      x0background[idx...] -= delta_x
      Ax0 = eval(p, x0background, param_local, idx, global_param)
      (Ax1 - Ax0)/(2delta_x)
    end

  end
  return Ai

end
export get_all_Ai

type Jaccobi_onestep{T <: Localized_Functor} <: GridMatrixIter_one
  A :: T
  OL :: Tuple{Vararg{Int64}}
end

function eval!{T<: Localized_Functor,n}(
  Ja :: Jaccobi_onestep{T}
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
  change_x = fill(0.0,size(x))

  Ais =  get_all_Ai(
           Ja.OL, Ja.A, x, param, global_param; 
           iseq = iseq, delta_x = delta_x, OL_bound = Ja.OL
           )
  updatepoint_range = ( 0 == length(iseq)) ? collect(1:reduce(*,irank)) : iseq
  for (ilogger, ipoint) = enumerate(updatepoint_range)
    iixy = ind2sub(irank,ipoint)
    iixyrange = map((1:n...)) do ii 
      OL[ii] + ((iixy[ii]-OL[ii]):(iixy[ii]+OL[ii]))
    end

    param_local = map(param) do xx
      xx[iixyrange...]
    end

    xbackground = x[iixyrange...]
    idx = map(x->x+1, OL)

    Ax = eval(Ja.A, xbackground, param_local, idx, global_param)
    Ai = Ais[ilogger]

    y_localerror = y[ipoint] - Ax

    xbackground_idx = map(ii -> OL[ii] + iixy[ii], (1:n...))
    change_x[xbackground_idx...] += alpha*y_localerror / Ai
  end
  x[:] =  x+change_x
  return x
end

export Jaccobi_onestep
