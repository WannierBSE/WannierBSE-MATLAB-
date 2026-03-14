function epslnq = epslion_Q2D(q_data, epsQ2D_data, q)
    epslnq = zeros(size(q)); % preallocate to match q's shape

    % Below three logical masks cover all cases for each q element
    t_smaller = q < min(q_data);
    t_middle  = q >= min(q_data) & q < max(q_data(end-1));
    t_larger  = q >= max(q_data(end-1));

    epslnq(t_smaller) = 1;
    epslnq(t_middle)  = interp1(q_data, epsQ2D_data, q(t_middle), 'spline');
    epslnq(t_larger)  = interp1(q_data, epsQ2D_data, q(t_larger), 'linear');
end
