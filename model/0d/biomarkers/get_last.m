function [last_AP_t, last_AP_var] = get_last(t, var, Tstart, Tend)

% t start
idx_start = 1;
while t(idx_start) < Tstart
    idx_start = idx_start + 1;
end

% idx_Ti end
idx_end = length(t);
while t(idx_end) > Tend
    idx_end = idx_end - 1;
end

last_AP_t   = t(idx_start:idx_end);
last_AP_var = var(idx_start:idx_end);

end