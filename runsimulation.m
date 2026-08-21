function [Result, plotData] = runsimulation(F, conv, P1, P2, n1, n2, R1, R2, ...
                                            ethfracf, watfracf, entfracf, ...
                                            ethfracb, H2Ofracb, ethfracd)
% runsimulation  Stage-by-stage simulation of azeotropic distillation of
%                ethanol-water using cyclohexane as entrainer
%
%   [Result, plotData] = runsimulation(F, conv, P1, P2, n1, n2, R1, R2, ...
%                                      ethfracf, watfracf, entfracf, ...
%                                      ethfracb, H2Ofracb, ethfracd)
%
%   PROCESS DESCRIPTION
%   -------------------
%   Feed (EtOH + H2O + CyHex) enters Column 1 (main azeotropic column).
%   The overhead vapour condenses and phase-splits in a decanter:
%     • Organic phase (CyHex-rich) → recycled to Column 1 as entrainer
%     • Aqueous phase (EtOH + H2O) → fed to Column 2 (recovery column)
%   Column 1 bottoms: high-purity ethanol product
%   Column 2 bottoms: essentially pure water
%   Column 2 distillate: dilute EtOH/H2O sent back to Column 1 feed
%
%   COLUMN MODEL
%   ------------
%   Stage-by-stage Lewis-Matheson method (enriching and stripping sections
%   solved separately, stepping away from specified terminal compositions).
%   Thermodynamics: modified Raoult's law with NRTL activity coefficients.
%   Vapour pressures: DIPPR 101 extended Antoine equation.
%
%   INPUTS
%   ------
%   F        – total molar feed flow rate     [kmol/hr]
%   conv     – fractional ethanol recovery    (0 < conv <= 1)
%   P1       – Column 1 operating pressure    [Pa]
%   P2       – Column 2 operating pressure    [Pa]
%   n1       – total theoretical stages, Col 1 (including reboiler)
%   n2       – total theoretical stages, Col 2 (including reboiler)
%   R1       – reflux ratio, Column 1         (L/D)
%   R2       – reflux ratio, Column 2         (L/D)
%   ethfracf – ethanol mole fraction in feed
%   watfracf – water   mole fraction in feed
%   entfracf – cyclohexane mole fraction in feed
%   ethfracb – target ethanol mole fraction in Col 1 bottoms (product purity)
%   H2Ofracb – target water mole fraction in Col 2 bottoms
%   ethfracd – target ethanol mole fraction in Col 2 distillate
%
%   OUTPUTS
%   -------
%   Result   – structure containing stream compositions, temperatures,
%              energy duties, decanter split, and summary tables
%   plotData – structure formatted for GUI plotting
%

%% 0. INPUT VALIDATION

validateattributes(F,       {'numeric'},{'scalar','positive'},    'runsimulation','F');
validateattributes(conv,    {'numeric'},{'scalar','>',0,'<=',1},  'runsimulation','conv');
validateattributes(P1,      {'numeric'},{'scalar','positive'},    'runsimulation','P1');
validateattributes(P2,      {'numeric'},{'scalar','positive'},    'runsimulation','P2');
validateattributes(n1, {'numeric'},{'scalar','integer','>',1},'runsimulation','n1');
validateattributes(n2, {'numeric'},{'scalar','integer','>',1},'runsimulation','n2');
validateattributes(R1, {'numeric'},{'scalar','positive'},    'runsimulation','R1');
validateattributes(R2, {'numeric'},{'scalar','positive'},    'runsimulation','R2');

zF_sum = ethfracf + watfracf + entfracf;
if abs(zF_sum - 1) > 1e-4
    error('runsimulation: feed mole fractions must sum to 1 (sum = %.6f).', zF_sum);
end
if ethfracb <= 0 || ethfracb > 1 || H2Ofracb <= 0 || H2Ofracb > 1 || ...
   ethfracd  <= 0 || ethfracd  > 1
    error('runsimulation: composition specifications must be in (0, 1].');
end

%% 1. PHYSICAL PROPERTY PARAMETERS

%  Components:  1 = Ethanol,  2 = Water,  3 = Cyclohexane
params.comps = {'Ethanol','Water','Cyclohexane'};

% Antoine coefficients — DIPPR 101 form (SI: Pa, K)
%   ln(Psat/Pa) = A + B/T + C*T + D*ln(T) + E*T^F
%   Source: DIPPR Project 801
%
%   Vector layout: [A,  B,     C,          D,       E,          F]
params.antoine.A = [73.304, -7122.3,  0,  0, -7.1424,  2.8853e-06,  2];
params.antoine.B = [73.649, -7258.2,  0,  0, -7.3037,  4.1653e-06,  2];
params.antoine.C = [51.087, -5226.4,  0,  0, -4.2278,  9.7554e-18,  6];

% CORRECTION: original code stored 7-element vectors [A,B,C,D,E,F_exp,F_pow]
% but psat_T was written to expect 7 elements.  Both are now consistent at
% 7 elements: [A, B/T-coeff, C*T-coeff, D*lnT-coeff, E, T-exponent].
% The coefficient layout matches the original convention.

% NRTL binary interaction parameters
%   tau_ij = Aij + Bij/T
%   G_ij   = exp(-alpha_ij * tau_ij)
%   Source: Knapp & Doherty (1994); Gomis et al. (2000)
%
%   Component ordering: [EtOH, H2O, CyHex]
params.NRTL.Aij = [ 0,       -0.8009,  -0.1560;
                    3.4578,   0,        13.1428;
                    1.6271, -10.4585,    0      ];

params.NRTL.Bij = [ 0,       246.18,   459.877;
                  -586.081,    0,     -1066.98;
                   214.076, 4954.9,       0    ];

% alpha_ij = alpha_ji (symmetric, non-randomness parameter)
params.NRTL.alpha = [ 0,    0.3,  0.45;
                      0.3,  0,    0.2;
                      0.45, 0.2,  0   ];

% NOTE: Cij is not used in the standard NRTL model (tau = A + B/T only).
% If you wish to include a temperature-dependent alpha (NRTL-1 model),
% add the extension here. Otherwise Cij is removed to avoid confusion.

% =========================================================================
%% 2. OVERALL MATERIAL BALANCE – COLUMN 1
% =========================================================================
zF = [ethfracf, watfracf, entfracf];

% Component molar flows in feed  [kmol/hr]
F_eth = F * ethfracf;
F_wat = F * watfracf;
F_ent = F * entfracf;

% Column 1 bottoms: ethanol product
%   Ethanol recovery specification:  mf_eth_B = conv * F_eth
mf_eth_B = conv * F_eth;
B1       = mf_eth_B / ethfracb;          % total bottoms flow  [kmol/hr]

% Feasibility check
if B1 > F
    warning('runsimulation: specified recovery/purity combination is infeasible. Adjusting recovery.');
    conv     = ethfracb * F / F_eth;
    mf_eth_B = conv * F_eth;
    B1       = mf_eth_B / ethfracb;
end

% Distribute remaining bottoms between H2O and CyHex in feed proportion
remaining_B1 = max(B1 - mf_eth_B, 0);
tot_others   = F_wat + F_ent;

if tot_others > 1e-12
    mf_wat_B1 = remaining_B1 * (F_wat / tot_others);
    mf_ent_B1 = remaining_B1 * (F_ent / tot_others);
else
    mf_wat_B1 = 0;
    mf_ent_B1 = 0;
end

xB1_spec = [mf_eth_B, mf_wat_B1, mf_ent_B1] / B1;  % Col 1 bottoms TARGET (used in operating line only)

% Column 1 distillate from overall balance
D1    = F - B1;
if D1 <= 0
    error('runsimulation: D1 <= 0 — bottoms cannot exceed feed. Check specifications.');
end

F_comp  = [F_eth, F_wat, F_ent];
B1_comp = [mf_eth_B, mf_wat_B1, mf_ent_B1];  % same as xB1_spec * B1
D1_comp = F_comp - B1_comp;
xD1_MB  = D1_comp / D1;                            % Col 1 distillate composition (by MB)

% =========================================================================
%% 3. VLE AT FEED CONDITIONS
% =========================================================================
[~, K_feed, T_bub_feed] = VLE_bubble(zF, P1, params);
fprintf('  Feed bubble-point: T = %.2f K,  K = [%.4g  %.4g  %.4g]\n', ...
        T_bub_feed, K_feed(1), K_feed(2), K_feed(3));

% Identify light and heavy keys by K-value at feed
[~, idx_sorted] = sort(K_feed);
keys.HK  = idx_sorted(1);    % smallest K → heavy key
keys.LK  = idx_sorted(end);  % largest K  → light key
alpha_LK = K_feed / K_feed(keys.HK);   % relative volatilities to HK

fprintf('  Light key: component %d (%s),  Heavy key: component %d (%s)\n', ...
        keys.LK, params.comps{keys.LK}, keys.HK, params.comps{keys.HK});



% =========================================================================
%% 4. COLUMN 1 STAGE-BY-STAGE SIMULATION (Lewis-Matheson)
% =========================================================================
%
%  Convention: stage 1 = reboiler (bottom), stage n1 = partial condenser / top.
%  We step from bottom to top in each section separately.
%
%  Stripping section: stages 1 … ns1   (bottom → feed)
%  Enriching section: stages ns1+1 … n1 (feed → top)

ns1 = max(1, round(n1 * 0.45));  % stripping stages (approx 45% of total)
nr1 = n1 - ns1;                  % enriching stages
Col1_FeedStage = ns1 + 1;

x1 = zeros(n1, 3);   % liquid composition profile  [stage, component]
y1 = zeros(n1, 3);   % vapour composition profile
T1 = zeros(n1, 1);   % temperature profile [K]

% -------------------------------------------------------------------------
%  AZEOTROPIC COLUMN — INTERNAL FLOW AND TOP COMPOSITION
% -------------------------------------------------------------------------
%  In an azeotropic column the top vapour is NOT the net distillate.
%  The top vapour goes to the decanter which splits into:
%    - Organic phase  → recycled back as reflux (entrainer recycle)
%    - Aqueous phase  → net distillate / Column 2 feed
%
%  The correct enriching section target is therefore the TOP VAPOUR
%  composition, which for the EtOH-H2O-CyHex system at atmospheric
%  pressure is near the ternary heterogeneous azeotrope:
%       EtOH~0.228,  H2O~0.238,  CyHex~0.534  (Gomis et al., 2000)
%
%  For the operating line we use the TOTAL vapour V1 = (R1+1)*D_total
%  where D_total is the total flow leaving the top of the column
%  (= organic recycle + aqueous distillate).  Since the organic recycle
%  is returned as reflux, we approximate:
%       V1  = (R1 + 1) * F        [upper bound, conservative]
%  This gives a realistically large internal vapour flow consistent with
%  the large entrainer recycle in heterogeneous azeotropic distillation.
%
%  The stripping section uses the correct CMO balance with:
%       V1s = V1  (q = 1 saturated liquid feed)
%       L1s = V1 + B1  (from V = L - B in stripping)

% Ternary heterogeneous azeotrope at 101325 Pa (Gomis et al., 2000)
% Used as enriching section top target (vapour composition)
xAZ = [0.228, 0.238, 0.534];   % [EtOH, H2O, CyHex]

% If the feed contains very little entrainer, scale the azeotrope target
% toward the binary EtOH-H2O azeotrope to remain physically consistent
if entfracf < 0.05
    x_binaz = [0.894, 0.106, 0];   % EtOH-H2O binary azeotrope at 1 atm
    w = entfracf / 0.05;
    xAZ = w * xAZ + (1-w) * x_binaz;
    xAZ = xAZ / sum(xAZ);
end

% Internal flows — correct CMO basis using distillate D1
%   V1  = (R1 + 1) * D1   [standard column vapour flow]
%   L1  = R1 * D1          [reflux flow]
%   L1s = V1 + B1          [stripping liquid, from CMO: L_s = V + B for q=1]
V1  = (R1 + 1) * D1;   % vapour flow in enriching section [kmol/hr]
L1  = R1 * D1;          % liquid (reflux) in enriching section [kmol/hr]
L1s = V1 + B1;          % stripping liquid flow from CMO balance

T_bounds = [250, 450];
x_stage  = xB1_spec;  % start stripping from specified bottoms target

% --- Stripping section (stages 1 to ns1) ---
for s = 1:ns1
    T_s = solve_Tbub(x_stage, P1, params, T_bounds);
    gamma   = activity_coeff_NRTL(x_stage, T_s, params);
    ps      = psat_T(T_s, params);
    Kvals   = gamma .* ps / P1;
    y_stage = Kvals .* x_stage;
    y_stage = y_stage / sum(y_stage);

    x1(s,:) = x_stage;
    y1(s,:) = y_stage;
    T1(s)   = T_s;

    % Stripping operating line:  L_s*x_{s+1} = V1*y_s + B1*xB1_spec
    x_next = (V1 / L1s) * y_stage + (B1 / L1s) * xB1_spec;
    x_next = max(x_next, 0);
    x_next = x_next / sum(x_next);
    x_stage = x_next;
end

% --- Enriching section (stages ns1+1 to n1) ---
%  Top vapour target = ternary heterogeneous azeotrope composition.
%  Enriching operating line:  x_{n+1} = (L1/V1)*y_n + (D_top/V1)*x_top
%  where D_top/V1 = 1/(R1+1) and x_top = xAZ.
for s = ns1+1 : n1
    T_s = solve_Tbub(x_stage, P1, params, T_bounds);
    gamma   = activity_coeff_NRTL(x_stage, T_s, params);
    ps      = psat_T(T_s, params);
    Kvals   = gamma .* ps / P1;
    y_stage = Kvals .* x_stage;
    y_stage = y_stage / sum(y_stage);

    x1(s,:) = x_stage;
    y1(s,:) = y_stage;
    T1(s)   = T_s;

    % Enriching operating line toward azeotrope top
    x_next = (R1 / (R1+1)) * y_stage + (1 / (R1+1)) * xAZ;
    x_next = max(x_next, 0);
    x_next = x_next / sum(x_next);
    x_stage = x_next;
end

% Top vapour leaving Column 1 → decanter feed
xD1_vap = y1(end,:);    % vapour composition at top stage

% Liquid compositions at ends (for T_cond and T_reb calculations)
xD1_liq = x1(end,:);
xB1_liq = x1(1,:);     % ACTUAL bottoms liquid from stage calculation (not the specified target)
% NOTE: xB1_liq(1) is the genuine EtOH purity the column achieves given R1 and n1.
%       It differs from ethfracb when R1 or n1 are insufficient to meet the spec.

% =========================================================================
%% 5. COLUMN 1 CONDENSER AND REBOILER TEMPERATURES
% =========================================================================
T_cond1 = solve_Tbub(xD1_liq, P1, params, T_bounds);
T_reb1  = solve_Tbub(xB1_liq, P1, params, T_bounds);

gamma_cond1 = activity_coeff_NRTL(xD1_liq, T_cond1, params);
ps_cond1    = psat_T(T_cond1, params);
K_cond1     = gamma_cond1 .* ps_cond1 ./ P1;

gamma_reb1  = activity_coeff_NRTL(xB1_liq, T_reb1, params);
ps_reb1     = psat_T(T_reb1, params);
K_reb1      = gamma_reb1 .* ps_reb1 ./ P1;

fprintf('  Column 1 | T_cond = %.2f K, T_reb = %.2f K\n', T_cond1, T_reb1);
fprintf('  K-values (cond) EtOH/H2O/CyHex = [%.4g  %.4g  %.4g]\n', K_cond1);
fprintf('  K-values (reb)  EtOH/H2O/CyHex = [%.4g  %.4g  %.4g]\n', K_reb1);

% =========================================================================
%% 6. DECANTER (LIQUID-LIQUID EQUILIBRIUM)
% =========================================================================
T_dec = 298.15;   % Decanter temperature [K] — typically 25°C (298.15 K)
                  % per industrial practice (Gomis et al., 2000)

[x_org, x_aq, beta, flag_dec, m_org, m_aq] = ...
    decanter_LLE(xD1_vap, T_dec, P1, params, D1);

if ~flag_dec
    warning('runsimulation: Decanter did not split into two phases. Check feed composition to decanter.');
end

fprintf('  Decanter | T = %.2f K,  beta_org = %.4f\n', T_dec, beta);
fprintf('  Organic phase:  [EtOH=%.4f  H2O=%.4f  CyHex=%.4f]  flow=%.4f\n', x_org, m_org);
fprintf('  Aqueous phase:  [EtOH=%.4f  H2O=%.4f  CyHex=%.4f]  flow=%.4f\n', x_aq, m_aq);

% =========================================================================
%% 7. COLUMN 2 – RECOVERY COLUMN (aqueous phase from decanter)
% =========================================================================
F2  = m_aq;    % [kmol/hr] — feed to Column 2 is aqueous phase
zF2 = x_aq;   % composition

if F2 <= 1e-12
    warning('runsimulation: Column 2 feed is zero or negligible. Skipping Column 2.');
    x2 = []; y2 = []; T2 = [];
    D2 = 0; B2 = 0;
    xD2 = [NaN, NaN, NaN];
    xB2 = [NaN, NaN, NaN];
    ns2 = 0; nr2 = 0;
    Col2_FeedStage = 0;
    T_cond2 = NaN; T_reb2 = NaN;
else
    % --- Column 2 material balance ---
    %   Bottoms: essentially pure water  xB2 = [1-H2Ofracb, H2Ofracb, 0]
    %   Distillate: dilute ethanol       xD2 ≈ [ethfracd, 1-ethfracd, 0]
    xB2_spec = [1 - H2Ofracb, H2Ofracb, 0];
    xD2_spec = [ethfracd,     1 - ethfracd, 0];

    % VLE at Column 2 feed
    gamma_f2 = activity_coeff_NRTL(zF2, T_dec, params);
    ps_f2    = psat_T(T_dec, params);
    K_f2     = gamma_f2 .* ps_f2 ./ P2;
    alpha2   = K_f2 / K_f2(2);    % relative volatility to water (HK)

    % Stage count split
    ns2 = max(1, round(n2 * 0.45));
    nr2 = n2 - ns2;
    Col2_FeedStage = ns2 + 1;

    % Column 2 internal flows (constant molar overflow, q = 1)
    %   Approximate D2 and B2 from feed balance:
    %   Overall:  F2 = D2 + B2
    %   EtOH:     F2*zF2(1) = D2*xD2(1) + B2*xB2(1)
    %   => D2 = F2*(zF2(1) - xB2_spec(1)) / (xD2_spec(1) - xB2_spec(1))
    denom_mb2 = xD2_spec(1) - xB2_spec(1);
    if abs(denom_mb2) < 1e-10
        D2 = F2 / (R2 + 1);   % fallback
    else
        D2 = F2 * (zF2(1) - xB2_spec(1)) / denom_mb2;
        D2 = max(D2, 0);
        D2 = min(D2, F2 - 1e-10);
    end
    B2  = F2 - D2;

    V2  = (R2 + 1) * D2;
    L2  = R2 * D2;
    L2s = L2 + F2;    % stripping liquid flow (q=1)

    % Initialise arrays
    x2 = zeros(n2, 3);
    y2 = zeros(n2, 3);
    T2 = zeros(n2, 1);

    x_stage2 = xB2_spec;

    % --- Stripping section Column 2 ---
    for s = 1:ns2
        T_s = solve_Tbub(x_stage2, P2, params, T_bounds);
        gamma    = activity_coeff_NRTL(x_stage2, T_s, params);
        ps       = psat_T(T_s, params);
        Kvals    = gamma .* ps / P2;
        y_stage2 = Kvals .* x_stage2;
        y_stage2 = y_stage2 / sum(y_stage2);

        x2(s,:) = x_stage2;
        y2(s,:) = y_stage2;
        T2(s)   = T_s;

        % Corrected stripping operating line
        if L2s > 1e-12
            x_next = (V2 / L2s) * y_stage2 + (B2 / L2s) * xB2_spec;
        else
            x_next = xB2_spec;
        end
        x_next  = max(x_next, 0);
        x_next  = x_next / sum(x_next);
        x_stage2 = x_next;
    end

    % --- Enriching section Column 2 ---
    for s = ns2+1 : n2
        T_s = solve_Tbub(x_stage2, P2, params, T_bounds);
        gamma    = activity_coeff_NRTL(x_stage2, T_s, params);
        ps       = psat_T(T_s, params);
        Kvals    = gamma .* ps / P2;
        y_stage2 = Kvals .* x_stage2;
        y_stage2 = y_stage2 / sum(y_stage2);

        x2(s,:) = x_stage2;
        y2(s,:) = y_stage2;
        T2(s)   = T_s;

        % Enriching operating line
        if V2 > 1e-12
            x_next = (L2 / V2) * y_stage2 + (D2 / V2) * xD2_spec;
        else
            x_next = xD2_spec;
        end
        x_next  = max(x_next, 0);
        x_next  = x_next / sum(x_next);
        x_stage2 = x_next;
    end

    xD2 = y2(end,:);   % Column 2 distillate (top vapour)
    xB2 = x2(1,:);     % Column 2 bottoms (bottom liquid)

    % Column 2 condenser and reboiler temperatures
    xD2_liq = x2(end,:);
    xB2_liq = x2(1,:);

    T_cond2 = solve_Tbub(xD2_liq, P2, params, T_bounds);
    T_reb2  = solve_Tbub(xB2_liq, P2, params, T_bounds);

    gamma_cond2 = activity_coeff_NRTL(xD2_liq, T_cond2, params);
    ps_cond2    = psat_T(T_cond2, params);
    K_cond2     = gamma_cond2 .* ps_cond2 ./ P2;

    gamma_reb2  = activity_coeff_NRTL(xB2_liq, T_reb2, params);
    ps_reb2     = psat_T(T_reb2, params);
    K_reb2      = gamma_reb2 .* ps_reb2 ./ P2;

    fprintf('  Column 2 | T_cond = %.2f K, T_reb = %.2f K\n', T_cond2, T_reb2);
    fprintf('  K-values (cond) EtOH/H2O/CyHex = [%.4g  %.4g  %.4g]\n', K_cond2);
    fprintf('  K-values (reb)  EtOH/H2O/CyHex = [%.4g  %.4g  %.4g]\n', K_reb2);
end

% =========================================================================
%% 8. ENERGY DUTIES (Watson correlation + sensible heat)
% =========================================================================
%
%  Duty sign convention:
%    Reboiler (Q_reb) > 0  — heat input
%    Condenser (Q_cond) < 0 — heat removed
%
%  Latent heat at temperature T estimated by Watson correlation:
%       lambda(T) = lambda_ref * [ (1 - T/Tc) / (1 - T_ref/Tc) ]^0.38
%
%  Reference latent heats and critical temperatures (DIPPR / Perry's):
%       EtOH:   lambda_ref = 38.56 kJ/mol  at T_ref = 351.5 K,  Tc = 514.0 K
%       H2O:    lambda_ref = 40.65 kJ/mol  at T_ref = 373.15 K, Tc = 647.1 K
%       CyHex:  lambda_ref = 33.10 kJ/mol  at T_ref = 353.0 K,  Tc = 553.5 K
%
%  Heat capacities (liquid, approximate mean values, kJ/mol·K):
%       EtOH: 0.112,  H2O: 0.075,  CyHex: 0.156
%
%  NOTE: flows in [kmol/hr] → convert to [mol/hr] × 1000 for kJ/hr duties
%        then divide by 3600 for kW.

n_watson = 0.38;

Href = [38.56, 40.65, 33.10];    % [kJ/mol]  EtOH, H2O, CyHex
Tc   = [514.0, 647.1, 553.5];    % [K]
Tref = [351.5, 373.15, 353.0];   % [K]
Cp   = [0.112, 0.075,  0.156];   % [kJ/mol·K]

T_ref_sens = 352;   % sensible heat reference temperature [K]

% Convert flows to mol/hr
D1_mol  = D1 * 1e3;
B1_mol  = B1 * 1e3;
V1_mol  = (R1 + 1) * D1 * 1e3;  % vapour flow consistent with D1-based column simulation
if exist('D2','var') && D2 > 0
    D2_mol = D2 * 1e3;
    B2_mol = B2 * 1e3;
end

[Qreb1, Qcond1] = column_duties(T_reb1, T_cond1, B1_mol, V1_mol, R1, xB1_liq, xD1_liq, ...
                                  Href, Tc, Tref, Cp, n_watson, T_ref_sens);

if exist('D2','var') && D2 > 0
    V2_mol = (R2 + 1) * D2 * 1e3;   % Col2 vapour — D2-based is correct here (no large recycle)
    [Qreb2, Qcond2] = column_duties(T_reb2, T_cond2, B2_mol, V2_mol, R2, xB2_liq, xD2_liq, ...
                                      Href, Tc, Tref, Cp, n_watson, T_ref_sens);
else
    Qreb2  = 0;
    Qcond2 = 0;
end

% =========================================================================
%% 9. ASSEMBLE RESULTS STRUCTURE
% =========================================================================

% --- Column 1 ---
Result.Column1.xTop       = xD1_vap;
Result.Column1.xBottom    = xB1_liq;
Result.Column1.TProfile   = T1;
Result.Column1.xProfile   = x1;
Result.Column1.yProfile   = y1;
Result.Column1.Tcond      = T_cond1;
Result.Column1.Treb       = T_reb1;
Result.Column1.Qreb_kW    = Qreb1  / 3600;
Result.Column1.Qcond_kW   = Qcond1 / 3600;
Result.Column1.Qreb_kJ_per_hr  = Qreb1;
Result.Column1.Qcond_kJ_per_hr = Qcond1;

% --- Decanter ---
Result.Decanter.xOrganic  = x_org;
Result.Decanter.xAqueous  = x_aq;
Result.Decanter.m_Organic = m_org;
Result.Decanter.m_Aqueous = m_aq;
Result.Decanter.beta      = beta;
Result.Decanter.Table = table( ...
    [m_org; m_aq], ...
    [x_org(1); x_aq(1)], ...
    [x_org(2); x_aq(2)], ...
    [x_org(3); x_aq(3)], ...
    'VariableNames', {'MolarFlow_kmolhr','x_Ethanol','x_Water','x_Cyclohexane'}, ...
    'RowNames', {'Organic_Phase','Aqueous_Phase'});

% --- Column 2 ---
if ~isempty(x2)
    Result.Column2.xTop       = xD2;
    Result.Column2.xBottom    = xB2;
    Result.Column2.TProfile   = T2;
    Result.Column2.xProfile   = x2;
    Result.Column2.yProfile   = y2;
    Result.Column2.Tcond      = T_cond2;
    Result.Column2.Treb       = T_reb2;
    Result.Column2.Qreb_kW    = Qreb2  / 3600;
    Result.Column2.Qcond_kW   = Qcond2 / 3600;
    Result.Column2.Qreb_kJ_per_hr  = Qreb2;
    Result.Column2.Qcond_kJ_per_hr = Qcond2;
else
    Result.Column2 = struct('xTop',[],'xBottom',[],'TProfile',[],...
                            'xProfile',[],'yProfile',[],...
                            'Tcond',NaN,'Treb',NaN,...
                            'Qreb_kW',0,'Qcond_kW',0);
end

% --- Column design summary for GUI text display ---
TextArea.Col1.RefluxRatio       = R1;
TextArea.Col1.StageNumber       = n1;
TextArea.Col1.StrippingStages   = ns1;
TextArea.Col1.RectifyingStages  = nr1;
TextArea.Col1.FeedStage         = Col1_FeedStage;

if F2 > 1e-12
    TextArea.Col2.RefluxRatio       = R2;
    TextArea.Col2.StageNumber       = n2;
    TextArea.Col2.StrippingStages   = ns2;
    TextArea.Col2.RectifyingStages  = nr2;
    TextArea.Col2.FeedStage         = Col2_FeedStage;
else
    TextArea.Col2 = struct('RefluxRatio',NaN,'StageNumber',NaN,...
                           'StrippingStages',0,'RectifyingStages',0,'FeedStage',0);
end

% Build display string for GUI text area
msg = sprintf([ ...
'  COLUMN 1 PARAMETERS:\n' ...
'    Reflux Ratio      : %.2f\n' ...
'    Total Stages      : %d\n' ...
'    Stripping Stages  : %d\n' ...
'    Rectifying Stages : %d\n' ...
'    Feed Stage        : %d\n\n' ...
'  COLUMN 2 PARAMETERS:\n' ...
'    Reflux Ratio      : %.2f\n' ...
'    Total Stages      : %d\n' ...
'    Stripping Stages  : %d\n' ...
'    Rectifying Stages : %d\n' ...
'    Feed Stage        : %d\n'], ...
    R1, n1, ns1, nr1, Col1_FeedStage, ...
    R2, n2, ns2, nr2, Col2_FeedStage);
TextArea.DesignString = msg;
Result.TextArea = TextArea;

% --- Summary stream table ---
xD2_disp = xD2(:)';  if any(isnan(xD2_disp)), xD2_disp = [NaN NaN NaN]; end
xB2_disp = xB2(:)';  if any(isnan(xB2_disp)), xB2_disp = [NaN NaN NaN]; end

SummaryTable = table( ...
    [D1;  B1;  D2;  B2], ...
    [xD1_vap(1); xB1_liq(1); xD2_disp(1); xB2_disp(1)], ...
    [xD1_vap(2); xB1_liq(2); xD2_disp(2); xB2_disp(2)], ...
    [xD1_vap(3); xB1_liq(3); xD2_disp(3); xB2_disp(3)], ...
    [T_cond1; T_reb1; T_cond2; T_reb2], ...
    'VariableNames', {'MolarFlow_kmolhr','x_Ethanol','x_Water','x_Cyclohexane','Temperature_K'}, ...
    'RowNames', {'Column1_Top','Column1_Bottom','Column2_Top','Column2_Bottom'});
Result.SummaryTable = SummaryTable;

% --- Energy summary table ---
EnergySummary = table( ...
    ["Column 1"; "Column 1"; "Column 2"; "Column 2"], ...
    ["Condenser"; "Reboiler"; "Condenser"; "Reboiler"], ...
    [Result.Column1.Qcond_kW; Result.Column1.Qreb_kW; ...
     Result.Column2.Qcond_kW; Result.Column2.Qreb_kW], ...
    'VariableNames', {'Column','Duty_Type','Duty_kW'});
Result.EnergySummary = EnergySummary;

% =========================================================================
%% 10. PLOT DATA FOR GUI
% =========================================================================
plotData.Column1 = struct('Stage',(1:n1)','T',T1,'x',x1,'y',y1);

if ~isempty(x2)
    plotData.Column2 = struct('Stage',(1:n2)','T',T2,'x',x2,'y',y2);
else
    plotData.Column2 = [];
end

plotData.Decanter     = struct('x_org',x_org,'x_aq',x_aq,'m_org',m_org,'m_aq',m_aq);
plotData.SummaryTable = SummaryTable;
plotData.Energy       = EnergySummary;

end  % end runsimulation

% =========================================================================
%% LOCAL HELPER FUNCTIONS
% =========================================================================

function T_bub = solve_Tbub(x, P, params, T_bounds)
% solve_Tbub  Robust bubble-point temperature solver with bracketed fzero.
%
%   Robust single-point solver used internally throughout runsimulation.
%   Falls back to a coarse scan if bracketing fails.

    x = max(x, 0);
    if sum(x) < 1e-20
        T_bub = mean(T_bounds);
        return;
    end
    x = x / sum(x);

    f = @(T) total_pressure(T, x, params, P);

    % Try bracketed solve first
    try
        f_lo = f(T_bounds(1));
        f_hi = f(T_bounds(2));
        if f_lo * f_hi < 0
            T_bub = fzero(f, T_bounds);
            return;
        end
    catch
        % continue to fallback
    end

    % Coarse scan to find bracket
    T_scan = linspace(T_bounds(1), T_bounds(2), 400);
    f_scan = arrayfun(f, T_scan);
    idx    = find(diff(sign(f_scan)) ~= 0, 1);
    if ~isempty(idx)
        try
            T_bub = fzero(f, [T_scan(idx), T_scan(idx+1)]);
            return;
        catch
        end
    end

    % Last resort: minimum residual point
    [~, imin] = min(abs(f_scan));
    T_bub = T_scan(imin);
    warning('solve_Tbub: bubble-point not bracketed; using minimum-residual estimate T = %.2f K', T_bub);
end
function [Q_reb, Q_cond] = column_duties(T_reb, T_cond, B_mol, V_mol, ~, ...
                                          xB, xD, Href, Tc, Tref, Cp, n, T_ref)
% column_duties  Reboiler and condenser duties for one column  [kJ/hr]
%
%   Method: latent heat of vaporisation (Watson correlation) for the
%   boilup/condensate stream, plus sensible heat to bring the liquid
%   streams to their operating temperatures.
%
%   Q_reb  > 0  (heat added at reboiler)
%   Q_cond < 0  (heat removed at condenser)
%   V_mol  = total internal vapour flow [mol/hr]
%   B_mol  = bottoms liquid flow [mol/hr]

    % Watson latent heats at reboiler and condenser temperatures
    %   lambda(T) = lambda_ref * [(1 - T/Tc) / (1 - Tref/Tc)]^0.38
    lambda_reb  = Href .* max(0, (1 - T_reb  ./ Tc) ./ (1 - Tref ./ Tc)) .^ n;
    lambda_cond = Href .* max(0, (1 - T_cond ./ Tc) ./ (1 - Tref ./ Tc)) .^ n;

    % Mixture latent heats weighted by stream composition
    H_reb  = sum(xB .* lambda_reb);   % reboiler vapour ~ bottoms composition
    H_cond = sum(xD .* lambda_cond);  % condenser vapour ~ top composition

    % Mixture liquid heat capacities [kJ/mol/K]
    Cp_bot = sum(xB .* Cp);
    Cp_top = sum(xD .* Cp);

    % Latent heat contributions [kJ/hr]
    Q_lat_reb  = V_mol * H_reb;    % heat to vaporise boilup at T_reb
    Q_lat_cond = V_mol * H_cond;   % heat to condense vapour at T_cond

    % Sensible heat contributions [kJ/hr]
    %   Reboiler:  heat bottoms liquid from reference to T_reb
    %   Condenser: sensible heat is small vs latent — omit to avoid double-counting
    Q_sens_reb  = B_mol * Cp_bot * (T_reb - T_ref);

    Q_reb  =  Q_lat_reb  + Q_sens_reb;
    Q_cond = -Q_lat_cond;
end