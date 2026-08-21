function dp = dpsat_dT(T, params)
% dpsat_dT  Analytical derivative of vapour pressure wrt temperature
%
%   dp = dpsat_dT(T, params)
%
%   Differentiates the DIPPR 101 expression analytically:
%       ln(Psat) = A + B/T + C*T + D*ln(T) + E*T^F
%       d(Psat)/dT = Psat * ( -B/T^2 + C + D/T + E*F*T^(F-1) )
%
%   This replaces the previous numerical central-difference approximation,
%   which used a step h = 0.1 K — appropriate but less efficient and
%   potentially inaccurate near phase boundaries.
%
%   Input:
%       T      – scalar temperature [K]
%       params – structure with field params.antoine.{A,B,C} (each 1x6)
%
%   Output:
%       dp     – 1x3 row vector  dPsat/dT  [Pa/K]

    if ~isscalar(T)
        error('dpsat_dT: T must be a scalar.');
    end

    psat_vals = psat_T(T, params);   % reuse validated psat_T
    dp        = zeros(1, 3);

    coeff_cell = {params.antoine.A, params.antoine.B, params.antoine.C};

    for ic = 1:3
        c  = coeff_cell{ic};   % [A, B, C, D, E, F_coeff, n_exp]  (7 elements)
        % d(lnPsat)/dT = -B/T^2 + C + D/T + F_coeff * n_exp * T^(n_exp - 1)
        dln_dT = -c(2)/T^2 + c(3) + c(5)/T + c(6)*c(7)*T^(c(7)-1);
        dp(ic) = psat_vals(ic) * dln_dT;
    end
end







% function dp = dpsat_dT(T, params_local)
%     % Numerical derivative of vapor pressure wrt T
%     % Input:  T [K] (scalar or vector)
%     %         params_local with Antoine coefficients
%     % Output: dp [Pa/K] (same size as T)
% 
%     h = 1e-1;   % finite difference step (K)
% 
%     psat_plus  = psat_T(T + h, params_local);
%     psat_minus = psat_T(T - h, params_local);
% 
%     dp = (psat_plus - psat_minus) ./ (2*h);
% end