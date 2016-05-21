immutable AdamsBashforth5 <: TimeIntergral_onestep 
  rf :: Tuple{Vararg{Array{Float64,1}}}

  function AdamsBashforth5()
    rf1 = [1.0]
    rf2 = [3/2;-1/2]
    rf3 = [23/12; -4/3; 5/12]
    rf4 = [55/24; -59/24; 37/24; -3/8]
    rf5 = [1901/720; -1387/360; 109/30; -637/360; 251/720]
    rf = (rf1, rf2,rf3,rf4,rf5) 
    new(rf)
  end
end

function eval!(tf :: AdamsBashforth5, p :: ModelParams,
  data :: Tuple, trends :: Tuple,
  t :: Float64, dt :: Float64)

  new_trend = fTrend(data, p,t)
  all_trends = (new_trend, trends...)
  ab_order = length(all_trends)
  for ii = 1:ab_order
    delta_ii = map(all_trends[ii]) do f
      f* (p.rf[ii] * dt)
    end
    data = map(+, data, delta_ii)
  end
  TimeForward!(p,t,dt)
  return (data,new_trend)
end

