function val = total_pressure_error(T, x, params_local, P_local)
    % Same logic as total_pressure, but matches your VLE_bubble signature
    gamma = activity_coeff_NRTL(x, T, params_local);
    ps    = psat_T(T, params_local);
    val = sum(x .* gamma .* ps) - P_local;
end