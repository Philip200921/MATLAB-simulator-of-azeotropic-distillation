function val = total_pressure(T, x, params, P)
% total_pressure  Residual of modified Raoult's law for bubble-point
%
%   val = total_pressure(T, x, params, P)
%
%   Computes:  P_calc(T,x) - P_system
%
%   where  P_calc = sum_i [ x_i * gamma_i(T,x) * Psat_i(T) ]
%
%   Used as the objective for fzero / fsolve to find bubble-point T.
%
%   Input:
%       T      – scalar temperature [K]
%       x      – 1x3 liquid mole fraction vector (must sum to 1)
%       params – parameter structure (.NRTL, .antoine)
%       P      – system pressure [Pa]
%
%   Output:
%       val    – scalar residual [Pa]  (= 0 at bubble point)

    % --- guard against non-physical compositions ----------------------------
    x = max(x, 0);
    x = x / sum(x);

    gamma  = activity_coeff_NRTL(x, T, params);   % 1x3
    ps     = psat_T(T, params);                    % 1x3  [Pa]
    P_calc = sum(x .* gamma .* ps);
    val    = P_calc - P;
end




% function val = total_pressure(T, x, params, P)
%     % Computes P_calc - P_system  (for root finding)
%     % T: temperature (K)
%     % x: liquid composition (row vector)
%     % params: structure with .NRTL and .antoine
%     % P: system pressure (Pa)
% 
%     gamma = activity_coeff_NRTL(x, T, params);
%     ps    = psat_T(T, params);
%     P_calc = sum(x .* gamma .* ps);
%     val = P_calc - P;
% end