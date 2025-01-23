function ts = get_time_variable(trange, dt1, dt2, twin, bcl)

ti = trange(1);  count = 1;
ts = nan(1,1e6);
ts(1) = trange(1);
while ti < trange(2)
    if mod(ti,bcl)<twin
        dt = dt1;
    else
        dt = dt2;
    end
    
    ti = round(ti + dt,4);
    ts(count+1) = ti;
    count = count + 1;
end
ts = ts(1:count);
end
