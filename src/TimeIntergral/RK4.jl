immutable RK4 <: TimeIntergral_onestep
end

immutable RK4_oneStage{ns} 
  n :: Int64
  function RK4_oneStage()
    @assert eltype(ns) <: Int64
    @assert 0<ns
    @assert ns<=4
    new(ns)
  end
end
export RK4_oneStage

#= function TimeForward!( p :: ModelParams, t :: Float64, dt :: Float64) =#
#=   error("No function for abstract types") =#
#= end =#


#= function fTrend( data:: Tuple, p:: ModelParams) =#
#=   error("No function for abstract types") =#
#= end =#

function eval(tf :: RK4_oneStage{1}, p :: ModelParams, y :: Tuple, y_0 :: Tuple, dt :: Float64)
  k1 = fTrend(y, p)
  y1 = map((1:length(y)...)) do ii
    k1[ii]*dt/2 + y_0[ii]
  end
  (k1,y1,dt/2)
end
function eval(tf :: RK4_oneStage{2}, p :: ModelParams, y :: Tuple, y_0 :: Tuple, dt :: Float64)
  k1 = fTrend(y, p)
  y1 = map((1:length(y)...)) do ii
    k1[ii]*dt/2 + y_0[ii]
  end
  (k1,y1,0.0)
end
function eval(tf :: RK4_oneStage{3}, p :: ModelParams, y :: Tuple, y_0 :: Tuple, dt :: Float64)
  k1 = fTrend(y, p)
  y1 = map((1:length(y)...)) do ii
    k1[ii]*dt + y_0[ii]
  end
  (k1,y1,dt/2)
end
function eval(tf :: RK4_oneStage{4}, p :: ModelParams, y :: Tuple, y_0 :: Tuple, dt :: Float64)
  k1 = fTrend(y, p)
  y1 = y_0
  return (k1,y1,0.0)
end

function eval!(tf :: RK4, p :: ModelParams,
  data :: Tuple, trends :: Tuple,
  t :: Float64, dt :: Float64)

  (k1,y1,forwardDt1) = eval(RK4_oneStage{1}(), p , data, data, dt)
  TimeForward!(p, t, forwardDt1)
  t = t + forwardDt1

  (k2,y2,forwardDt2) = eval(RK4_oneStage{2}(), p , y1, data, dt)
  TimeForward!(p, t, forwardDt2)
  t = t + forwardDt2

  (k3,y3,forwardDt3) = eval(RK4_oneStage{3}(), p , y2, data, dt)
  TimeForward!(p, t, forwardDt3)
  t = t + forwardDt3

  (k4,y4,forwardDt4) = eval(RK4_oneStage{4}(), p , y3, data, dt)
  #= TimeForward!(p, t, forwardDt4) =#
  #= t = t + forwardDt4 =#

  yfinal = map((1:length(y1)...)) do ii
    data[ii] + (k1[ii]+2*k2[ii] + 2*k3[ii] +k4[ii])*dt/6
  end
  return (p, yfinal, t)
end
export RK4
