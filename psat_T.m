function psat = psat_T(T, params)
% psat_T  Compute pure-component vapour pressures via extended Antoine equation
%
%   psat = psat_T(T, params)
%
%   Model (DIPPR 101 form, units: Pa, K):
%       ln(Psat) = A + B/T + C*T + D*ln(T) + E*T^F
%
%   Components:  1 = Ethanol,  2 = Water,  3 = Cyclohexane
%
%   Parameters sourced from DIPPR 801 database (SI units: Pa, K).
%   Antoine coefficient vectors stored in params.antoine.A/B/C must have
%   exactly 7 elements:
%       [A, B, C, D, E, F, n]
%   where the power term is  E * T^n  (i.e., element 6 = E coefficient,
%   element 7 = n exponent).  This matches the original parameter set in
%   runsimulation.m:
%       Ethanol:    [73.304, -7122.3, 0, 0, -7.1424, 2.8853e-06, 2]
%       Water:      [73.649, -7258.2, 0, 0, -7.3037, 4.1653e-06, 2]
%       Cyclohexane:[51.087, -5226.4, 0, 0, -4.2278, 9.7554e-18, 6]
%
%   Input:
%       T      – scalar temperature [K]
%       params – structure with field params.antoine.{A,B,C} (each 1x7)
%
%   Output:
%       psat   – 1x3 row vector of saturation pressures [Pa]

    % --- input validation ---------------------------------------------------
    if ~isscalar(T)
        error('psat_T: T must be a scalar.');
    end
    if T < 200 || T > 600
        warning('psat_T: T = %.2f K is outside the expected range [200, 600] K.', T);
    end

    % --- evaluate for each component ----------------------------------------
    psat       = zeros(1, 3);
    coeff_cell = {params.antoine.A, params.antoine.B, params.antoine.C};

    for ic = 1:3
        c = coeff_cell{ic};

        if numel(c) ~= 7
            error('psat_T: Antoine coefficient vector for component %d must have 7 elements [A,B,C,D,E,F,n].', ic);
        end

        % ln(Psat) = A + B/T + C*T + D*ln(T) + E*T^n
        % c = [A, B, C, D, E, F_coeff, n_exp]
        %      1   2  3  4  5  6        7
        lnP      = c(1) + c(2)/T + c(3)*T + c(5)*log(T) + c(6)*(T^c(7));
        psat(ic) = exp(lnP);
    end
end



% function psat = psat_T(T, params)
% 
%     T = T(:);  % to ensure column vector
% 
%     % Component 1
%     A = params.antoine.A;
%     lnPsat1 = A(1) + A(2)./(T + A(3)) + A(4).*T + A(5).*log(T) + A(6).*(T.^A(7));
%     Psat1 = exp(lnPsat1);
% 
%     % Component 2
%     B = params.antoine.B;
%     lnPsat2 = B(1) + B(2)./(T + B(3)) + B(4).*T + B(5).*log(T) + B(6).*(T.^B(7));
%     Psat2 = exp(lnPsat2);
% 
%     % Component 3
%     C = params.antoine.C;
%     lnPsat3 = C(1) + C(2)./(T + C(3)) + C(4).*T + C(5).*log(T) + C(6).*(T.^C(7));
%     Psat3 = exp(lnPsat3);
% 
%     % Combine
%     psat = [Psat1, Psat2, Psat3];
% end
