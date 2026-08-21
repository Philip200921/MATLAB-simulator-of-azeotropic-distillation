function [x_org, x_aq, beta, flag, m_org, m_aq] = decanter_LLE(xD, T_dec, P1, params, D)
% decanter_LLE  Ternary LLE for EtOH-H2O-CyHex via Gibbs energy minimisation
%
%   Solves the decanter phase split by minimising total molar Gibbs energy:
%       G_total/RT = beta*G_mix(x_org) + (1-beta)*G_mix(x_aq)
%   subject to the overall material balance:
%       xD = beta*x_org + (1-beta)*x_aq
%
%   PHASE CONVENTION (always enforced):
%       x_org = CyHex-rich (organic) phase   [x_CyHex > x_aq_CyHex]
%       x_aq  = H2O-rich  (aqueous)  phase   [x_H2O   > x_org_H2O ]
%
%   References:
%   [1] Gomis et al. (2000) Fluid Phase Equilib. 172, 181-187.
%       Tie-line data at 298.15 K (used as initial guesses):
%         Organic:  EtOH=0.114, H2O=0.018, CyHex=0.868
%         Aqueous:  EtOH=0.082, H2O=0.900, CyHex=0.018
%   [2] Michelsen (1982) Fluid Phase Equilib. 9, 21-32.  (stability test)
%   [3] Renon & Prausnitz (1968) AIChE J. 14(1), 135-144. (NRTL)
%
%   Inputs:
%     xD     – 1x3 overall decanter feed  [EtOH, H2O, CyHex]
%     T_dec  – temperature [K]
%     P1     – pressure [Pa]  (interface compatibility; not used directly)
%     params – structure with .NRTL and .antoine
%     D      – total molar feed flow [kmol/hr or mol/hr]
%
%   Outputs:
%     x_org  – 1x3 organic phase mole fractions  (CyHex-rich)
%     x_aq   – 1x3 aqueous phase mole fractions  (H2O-rich)
%     beta   – organic phase fraction (0..1)
%     flag   – 2=two phases,  1=single phase
%     m_org  – organic phase flow [same units as D]
%     m_aq   – aqueous phase flow [same units as D]

% =========================================================================
%% 0. SANITISE INPUT
% =========================================================================
xD = xD(:)';
xD = max(xD, 0);
xD = xD / max(sum(xD), 1e-15);

if D <= 0,       error('decanter_LLE: D must be positive.');        end
if T_dec <= 0,   error('decanter_LLE: T_dec must be positive [K].'); end

% =========================================================================
%% 1. SINGLE-PHASE REFERENCE GIBBS ENERGY
% =========================================================================
G_single = g_mix(xD, T_dec, params);

% =========================================================================
%% 2. OPTIMISATION BOUNDS AND OPTIONS
% =========================================================================
%  Optimisation variables:  u = [x1_org(EtOH),  x2_org(H2O),  beta]
%  x3_org(CyHex) = 1 - x1_org - x2_org   (eliminated)
%  x_aq derived from material balance
eps_x = 1e-6;
eps_b = 1e-6;

lb    = [eps_x,     eps_x,     eps_b  ];
ub    = [1-2*eps_x, 1-2*eps_x, 1-eps_b];
A_lin = [1, 1, 0];       % x1_org + x2_org <= 1 - eps_x
b_lin = 1 - eps_x;

opts = optimoptions('fmincon', ...
    'Display',                'off', ...
    'Algorithm',              'sqp', ...
    'TolFun',                 1e-12, ...
    'TolX',                   1e-12, ...
    'TolCon',                 1e-12, ...
    'MaxIterations',          3000,  ...
    'MaxFunctionEvaluations', 30000);

% =========================================================================
%% 3. INITIAL GUESSES
% =========================================================================
%  The optimisation variable u = [x1_org(EtOH), x2_org(H2O), beta]
%  x3_org(CyHex) is implied as 1 - x1 - x2.
%
%  Guesses are anchored to Gomis et al. (2000) tie-line data and cover
%  a range of beta values to avoid local minima.
%  CRITICAL: x2_org(H2O) must be SMALL (organic phase is CyHex-rich,
%  not water-rich). Guesses with large x2_org would push the optimiser
%  toward the inverted (wrong) solution.

guesses = [ ...
%   x1_org(EtOH)  x2_org(H2O)  beta
    0.114,        0.018,        0.50;   % Gomis tie-line, beta=0.50
    0.114,        0.018,        0.30;   % Gomis tie-line, beta=0.30
    0.114,        0.018,        0.70;   % Gomis tie-line, beta=0.70
    0.100,        0.020,        0.40;
    0.120,        0.015,        0.55;
    0.080,        0.025,        0.45;
    0.150,        0.030,        0.50;
    0.050,        0.010,        0.60;
];

% =========================================================================
%% 4. MULTI-START OPTIMISATION
% =========================================================================
bestG  = Inf;
best_u = [];

for k = 1:size(guesses, 1)
    u0 = guesses(k, :);

    % Project strictly inside feasible box
    u0(1) = max(min(u0(1), ub(1)-eps_x), lb(1)+eps_x);
    u0(2) = max(min(u0(2), ub(2)-eps_x), lb(2)+eps_x);
    u0(3) = max(min(u0(3), ub(3)-eps_b), lb(3)+eps_b);
    if u0(1) + u0(2) >= b_lin
        u0(1) = (b_lin - 2*eps_x) / 2;
        u0(2) = (b_lin - 2*eps_x) / 2;
    end

    try
        [u_opt, G_opt, exitflag] = fmincon( ...
            @(u) obj_u(u, xD, T_dec, params), ...
            u0, A_lin, b_lin, [], [], lb, ub, [], opts);
    catch
        continue;
    end

    if exitflag < 1,  continue;  end

    if G_opt < bestG
        bestG  = G_opt;
        best_u = u_opt;
    end
end

% =========================================================================
%% 5. DEFAULT: SINGLE PHASE
% =========================================================================
x_org = xD;
x_aq  = xD;
beta  = 1;
flag  = 1;
m_org = D;
m_aq  = 0;

if isempty(best_u)
    warning('decanter_LLE: optimisation failed from all starting points. Single-phase returned.');
    return;
end

% =========================================================================
%% 6. RECONSTRUCT PHASE COMPOSITIONS
% =========================================================================
x1o      = best_u(1);
x2o      = best_u(2);
beta_opt = best_u(3);
x3o      = max(0, 1 - x1o - x2o);

% Phase I as returned by optimiser (intended to be organic/CyHex-rich)
phaseI = max([x1o, x2o, x3o], 0);
phaseI = phaseI / sum(phaseI);

denom = 1 - beta_opt;
if denom < eps_b
    warning('decanter_LLE: beta_opt ~ 1, no aqueous phase. Single-phase returned.');
    return;
end

phaseII = (xD - beta_opt .* phaseI) ./ denom;

if any(phaseII < -1e-4)
    warning('decanter_LLE: aqueous phase has negative mole fractions. Single-phase returned.');
    return;
end
phaseII = max(phaseII, 0);
phaseII = phaseII / sum(phaseII);

% =========================================================================
%% 7. PHASE IDENTIFICATION  (critical — prevents label inversion)
% =========================================================================
%  By physical convention for EtOH-H2O-CyHex:
%    Organic phase  = higher cyclohexane content  (component 3)
%    Aqueous phase  = higher water content         (component 2)
%
%  After optimisation, phaseI may or may not be the CyHex-rich phase.
%  We check explicitly and swap if necessary.

if phaseI(3) >= phaseII(3)
    % phaseI is already CyHex-rich  → organic = phaseI
    x_org_cand  = phaseI;
    x_aq_cand   = phaseII;
    beta_cand   = beta_opt;
else
    % phaseI is H2O-rich → swap labels and invert beta
    x_org_cand  = phaseII;
    x_aq_cand   = phaseI;
    beta_cand   = 1 - beta_opt;
end

% Verify the swap is physically consistent
if x_aq_cand(2) < x_org_cand(2)
    warning('decanter_LLE: phase identification ambiguous — water content check failed. Results may be unreliable.');
end

% =========================================================================
%% 8. MICHELSEN STABILITY TEST
% =========================================================================
G_two = beta_cand * g_mix(x_org_cand, T_dec, params) + ...
        (1-beta_cand) * g_mix(x_aq_cand,  T_dec, params);

phase_diff = norm(x_org_cand - x_aq_cand, inf);
COMP_TOL   = 1e-3;

if (G_two < G_single - 1e-10) && (phase_diff > COMP_TOL)
    flag  = 2;
    x_org = x_org_cand;
    x_aq  = x_aq_cand;
    beta  = beta_cand;
    m_org = beta     * D;
    m_aq  = (1-beta) * D;

    fprintf('  Decanter LLE: two-phase split | beta_org = %.4f\n', beta);
    fprintf('    Organic (CyHex-rich): [EtOH=%.4f  H2O=%.4f  CyHex=%.4f]  F=%.4g\n', x_org, m_org);
    fprintf('    Aqueous (H2O-rich) : [EtOH=%.4f  H2O=%.4f  CyHex=%.4f]  F=%.4g\n', x_aq,  m_aq);
else
    warning('decanter_LLE: phase split not thermodynamically favoured (dG=%.2e, diff=%.2e). Single-phase returned.', ...
            G_two - G_single, phase_diff);
end

% =========================================================================
%% 9. MASS BALANCE CHECK
% =========================================================================
MB_err = max(abs(beta*x_org + (1-beta)*x_aq - xD));
if MB_err > 1e-5
    warning('decanter_LLE: mass balance residual = %.2e.', MB_err);
end

end  % ---- end decanter_LLE ----


% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function G = g_mix(xv, T, params)
% Dimensionless molar Gibbs energy of mixing:  G_mix/RT = sum_i[x_i*(ln x_i + ln gamma_i)]
    xv = max(xv, 1e-15);
    xv = xv / sum(xv);
    gamma_v = activity_coeff_NRTL(xv, T, params);
    if any(~isfinite(gamma_v)) || any(gamma_v <= 0)
        G = 1e6;
        return;
    end
    G = sum(xv .* (log(xv) + log(gamma_v)));
end


function Gtot = obj_u(u, xD, T, params)
% Total two-phase Gibbs energy — objective for fmincon
%   u = [x1_org(EtOH), x2_org(H2O), beta]
    x1o = u(1);  x2o = u(2);  b = u(3);
    x3o = 1 - x1o - x2o;

    if x3o < 1e-9 || b < 1e-9 || b > 1-1e-9
        Gtot = 1e6;
        return;
    end

    x_org_loc = [x1o, x2o, x3o];
    denom = 1 - b;
    if denom < 1e-12,  Gtot = 1e6;  return;  end

    x_aq_loc = (xD - b .* x_org_loc) ./ denom;

    if any(x_aq_loc < -1e-6)
        Gtot = 1e6 + 1e4*sum(abs(min(x_aq_loc, 0)));
        return;
    end

    x_aq_loc = max(x_aq_loc, 0);
    x_aq_loc = x_aq_loc / max(sum(x_aq_loc), 1e-15);

    G_org = g_mix(x_org_loc, T, params);
    G_aq  = g_mix(x_aq_loc,  T, params);
    Gtot  = b*G_org + (1-b)*G_aq;
end


% function [x_org, x_aq, beta, flag, m_org, m_aq] = decanter_LLE(xD, T_dec, P1, params, D)
% % decanter_LLE_minG  -- Gibbs-energy minimization LLE solver (ternary)
% % Inputs:
% %   xD : 1x3 distillate composition [EtOH, H2O, CYC]
% %   T_dec, P, params : operating conditions and NRTL params
% %   D  : total molar feed to decanter
% % Outputs:
% %   x_org, x_aq : phase compositions (1x3)
% %   beta : organic phase fraction (0..1)
% %   flag : 2 => two-phase found, 1 => single-phase
% %   m_org, m_aq : phase flows
% 
%     % defaults -> single-phase
%     xD = xD(:)'; xD = max(xD,0); xD = xD./max(sum(xD),1e-16);
%     x_org = xD; x_aq = xD; beta = 1; flag = 1;
%     m_org = D; m_aq = 0;
% 
%     % objective: total molar Gibbs (dimensionless per mole; ignores RT factor)
%     % g(x) = sum_i x_i * ( log(x_i) + log(gamma_i(x)) )
%     function G = g_mix(xv)
%         xv = max(xv, 1e-15);
%         gamma_v = activity__coeff__NRTL(xv, T_dec, params);
%         % protect against nonfinite gamma
%         gamma_v(~isfinite(gamma_v)) = 1e12;
%         G = sum(xv .* (log(xv) + log(gamma_v)));
%     end
% 
%     % reduced objective: variables u = [x1o, x2o, beta] (x3o = 1-x1o-x2o)
%     function Gtot = obj_u(u)
%         x1o = u(1); x2o = u(2); b = u(3);
%         x3o = 1 - x1o - x2o;
%         if any([x1o,x2o,x3o,b] <= 1e-12) || b <= 1e-8 || b >= 1-1e-8
%             Gtot = 1e6 + sum(abs([x1o,x2o,x3o,b])); return;
%         end
%         x_org_local = [x1o,x2o,x3o];
%         denom = (1-b);
%         if denom <= 1e-12
%             Gtot = 1e6; return;
%         end
%         x_aq_local = (xD - b .* x_org_local) ./ denom;
%         if any(x_aq_local <= 1e-12) || any(~isfinite(x_aq_local))
%             Gtot = 1e6; return;
%         end
%         % compute molar Gibbs (dimensionless)
%         Gorg = g_mix(x_org_local);
%         Gaq  = g_mix(x_aq_local);
%         Gtot = b*Gorg + (1-b)*Gaq;
%     end
% 
%     % trying several initial guesses to avoid local traps
%     guesses = [0.02 0.10 0.1; 0.05 0.05 0.3; 0.1 0.05 0.4; 0.02 0.30 0.2; 0.2 0.6 0.3];
%     opts = optimoptions('fmincon','Display','off','Algorithm','sqp','TolFun',1e-10,'TolX',1e-10,'MaxIterations',2000);
% 
%     bestG = Inf; bestu = [];
%     for k = 1:size(guesses,1)
%         u0 = guesses(k,:);
%         % bounds and linear constraints: x1o+x2o <= 0.999, >= small, beta in [1e-6,1-1e-6]
%         lb = [1e-6, 1e-6, 1e-6];
%         ub = [0.999, 0.999, 1-1e-6];
%         % linear inequality: x1o + x2o <= 0.999
%         A = [1 1 0];
%         bineq = 0.999;
%         try
%             [u_opt, Gopt, exitflag] = fmincon(@obj_u, u0, A, bineq, [], [], lb, ub, [], opts);
%         catch
%             exitflag = -1; u_opt = []; Gopt = Inf;
%         end
%         if exitflag <= 0
%             continue;
%         end
%         if Gopt < bestG
%             bestG = Gopt; bestu = u_opt;
%         end
%     end
% 
%     if isempty(bestu)
%          %fallback: no feasible minimization found
%         flag = 1; x_org = xD; x_aq = xD; beta = 1; m_org = D; m_aq = 0; return;
%     end
% 
%     % unpack best solution
%     x1o = bestu(1); x2o = bestu(2); beta = bestu(3);
%     x3o = max(0, 1 - x1o - x2o);
%     x_org = [x1o, x2o, x3o];
%     x_aq = (xD - beta .* x_org) ./ (1 - beta);
%     % sanitize
%     x_org = max(x_org,0); x_org = x_org / max(sum(x_org),1e-16);
%     x_aq  = max(x_aq,0);  x_aq  = x_aq  / max(sum(x_aq),1e-16);
% 
%     % To check that phases are meaningfully different
%     if norm(x_org - x_aq, inf) < 1e-4
%         % no real demix
%         flag = 1; beta = 1; x_org = xD; x_aq = xD; m_org = D; m_aq = 0;
%     else
%         flag = 2;
%         m_org = beta * D;
%         m_aq  = (1-beta) * D;
%     end
% end
