
function [y, K, T_bub] = VLE_bubble(x, P, params)
% VLE_bubble  Bubble-point temperature and equilibrium vapour composition
%
%   [y, K, T_bub] = VLE_bubble(x, P, params)
%
%   Solves the bubble-point equation via modified Raoult's law (NRTL activity
%   coefficients):
%
%       sum_i [ x_i * gamma_i(T,x) * Psat_i(T) ] = P       ... (1)
%
%   Equation (1) is solved for T using fzero with a physically motivated
%   initial guess (mole-fraction-weighted boiling points).
%
%   Vapour compositions follow from:
%       y_i = x_i * gamma_i * Psat_i / P                    ... (2)
%   normalised to satisfy  sum(y_i) = 1.
%
%   Input:
%       x      – 1x3 liquid mole fraction row vector
%       P      – system pressure [Pa]
%       params – parameter structure (.NRTL, .antoine)
%
%   Output:
%       y      – 1x3 vapour mole fractions at bubble point
%       K      – 1x3 K-values  (K_i = y_i / x_i)
%       T_bub  – bubble-point temperature [K]

    % --- input checks -------------------------------------------------------
    x = x(:)';
    if numel(x) ~= 3
        error('VLE_bubble: only 3-component systems are supported.');
    end
    x = max(x, 0);
    x = x / sum(x);

    if P <= 0
        error('VLE_bubble: pressure P must be positive [Pa].');
    end

    % --- initial temperature guess (mole-weighted normal boiling points) ----
    % Ethanol 351.4 K, Water 373.1 K, Cyclohexane 353.9 K  (at 101325 Pa)
    Tb_pure = [351.4, 373.1, 353.9];
    T0      = max(x .* Tb_pure) * 3;   % weighted mean
    T0      = sum(x .* Tb_pure);

    % --- solve bubble-point equation ----------------------------------------
    f_bub = @(T) total_pressure(T, x, params, P);

    opts  = optimset('Display','off','TolFun',1e-9,'TolX',1e-9,'MaxIter',500);

    try
        T_bub = fzero(f_bub, [280, 450]);   % bracketed solve – more robust
    catch
        try
            T_bub = fzero(f_bub, T0, opts); % fallback: unbracketed from guess
        catch ME
            warning('VLE_bubble: fzero failed (%s). Using weighted boiling-point estimate.', ME.message);
            T_bub = T0;
        end
    end

    % --- compute K-values and vapour composition ----------------------------
    gamma  = activity_coeff_NRTL(x, T_bub, params);
    ps     = psat_T(T_bub, params);
    K      = gamma .* ps ./ P;

    Ky = K .* x;
    if sum(Ky) < 1e-20
        error('VLE_bubble: K*x sum is zero — check NRTL parameters and pressure.');
    end
    y = Ky / sum(Ky);    % normalise to unit sum
end

% end
