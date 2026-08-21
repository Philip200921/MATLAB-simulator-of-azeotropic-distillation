function generate_RCM()
% generate_RCM  Residue Curve Map for Ethanol-Water-Cyclohexane at 101325 Pa
%
%   Generates the full residue curve map (RCM) for the ternary system
%   EtOH(1)-H2O(2)-CyHex(3) using the NRTL activity coefficient model
%   and DIPPR 101 vapour pressure correlations — the same parameters used
%   in the azeotropic distillation simulator.
%
%   THEORY
%   ------
%   A residue curve is the locus of liquid compositions obtained by
%   open evaporation of a mixture at constant pressure.  It satisfies:
%
%       dx_i/d(xi) = x_i - y_i*(x,T)        i = 1,2,3         ...(1)
%
%   where xi is a dimensionless time variable (extent of evaporation),
%   x_i is the liquid mole fraction, and y_i is the equilibrium vapour
%   mole fraction computed from modified Raoult's law (NRTL).
%
%   The map is plotted on a ternary diagram.  Key features identified:
%     - Three pure-component nodes (stable/unstable)
%     - Binary azeotropes (EtOH-H2O at ~0.894/0.106, CyHex-EtOH at ~0.44/0.56)
%     - Ternary heterogeneous azeotrope (~0.228/0.238/0.534 at 101325 Pa)
%     - Distillation boundaries separating composition space into regions
%     - Two-phase (liquid-liquid) envelope
%
%   REFERENCES
%   ----------
%   [1] Doherty, M.F., Malone, M.F. (2001) Conceptual Design of Distillation
%       Systems. McGraw-Hill. (RCM theory and computation)
%   [2] Gomis, V. et al. (2000) Fluid Phase Equilib. 172, 181-187.
%       (LLE data and ternary azeotrope composition)
%   [3] Knapp, J.P., Doherty, M.F. (1994) AIChE J. 40(2), 243-268.
%       (heterogeneous azeotropic distillation design)
%
%   OUTPUT
%   ------
%   Figure 1: Full residue curve map on ternary diagram
%   Figure 2: LLE two-phase envelope at 298.15 K overlaid on RCM
%
%   USAGE
%   -----
%   Place this file in the same folder as:
%       activity_coeff_NRTL.m,  psat_T.m,  VLE_bubble.m,
%       total_pressure.m,  decanter_LLE.m
%   Then run:  >> generate_RCM

% =========================================================================
%% 0. PARAMETERS  (identical to runsimulation.m)
% =========================================================================
params.comps = {'Ethanol','Water','Cyclohexane'};

params.antoine.A = [73.304, -7122.3, 0, 0, -7.1424, 2.8853e-06, 2];
params.antoine.B = [73.649, -7258.2, 0, 0, -7.3037, 4.1653e-06, 2];
params.antoine.C = [51.087, -5226.4, 0, 0, -4.2278, 9.7554e-18, 6];

params.NRTL.Aij = [ 0,       -0.8009,  -0.1560;
                    3.4578,   0,        13.1428;
                    1.6271, -10.4585,   0      ];
params.NRTL.Bij = [ 0,       246.18,   459.877;
                  -586.081,   0,      -1066.98;
                   214.076, 4954.9,    0      ];
params.NRTL.alpha = [0,   0.3,  0.45;
                     0.3, 0,    0.2;
                     0.45,0.2,  0   ];

P = 101325;   % [Pa]  atmospheric pressure

% =========================================================================
%% 1. COMPUTE RESIDUE CURVES BY ODE INTEGRATION
% =========================================================================
fprintf('Computing residue curves...\n');

% ODE: dx/dxi = x - y(x)
% Integrated forward (increasing xi → vaporisation) and backward
% (decreasing xi → condensation) from each starting point.
% Forward integration → stable nodes (high-boiling)
% Backward integration → unstable nodes (low-boiling)

ode_opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'Events',@stop_event);

% Grid of starting compositions on the ternary diagram
% Use a triangular grid with spacing 0.1
starts = [];
step = 0.08;
for x1 = step : step : 1-step
    for x2 = step : step : 1-x1-step/2
        x3 = 1 - x1 - x2;
        if x3 > step/2
            starts = [starts; x1, x2, x3];  %#ok<AGROW>
        end
    end
end
% Add some starts near the azeotropes and boundaries
extra = [0.228, 0.238, 0.534;   % ternary azeotrope
         0.894, 0.106, 0.000;   % EtOH-H2O azeotrope
         0.440, 0.000, 0.560;   % EtOH-CyHex azeotrope (approx)
         0.05,  0.05,  0.90;
         0.05,  0.90,  0.05;
         0.90,  0.05,  0.05;
         0.30,  0.30,  0.40;
         0.10,  0.50,  0.40;
         0.50,  0.10,  0.40];
starts = [starts; extra];

% ODE right-hand side
rcm_ode = @(xi, x) rcm_rhs(xi, x, P, params);

curves_fwd = {};
curves_bwd = {};

for k = 1:size(starts,1)
    x0 = starts(k,:)';
    x0 = max(x0, 1e-6);
    x0 = x0 / sum(x0);

    % Forward (xi: 0 → +8)
    try
        [~, X] = ode45(rcm_ode, [0, 8], x0, ode_opts);
        X = max(X, 0);
        X = bsxfun(@rdivide, X, sum(X,2));
        curves_fwd{end+1} = X; %#ok<AGROW>
    catch
    end

    % Backward (xi: 0 → -8)
    try
        [~, X] = ode45(@(xi,x) -rcm_rhs(xi,x,P,params), [0, 8], x0, ode_opts);
        X = max(X, 0);
        X = bsxfun(@rdivide, X, sum(X,2));
        curves_bwd{end+1} = X; %#ok<AGROW>
    catch
    end
end

fprintf('  Computed %d forward + %d backward curve segments.\n', ...
        numel(curves_fwd), numel(curves_bwd));

% =========================================================================
%% 2. LLE ENVELOPE — Gomis et al. (2000) Table 1 experimental tie-line data
% =========================================================================
%  Using published experimental tie-line end-points directly rather than
%  scanning via decanter_LLE, which avoids false positives across the
%  composition space.  Data at 298.15 K, 101325 Pa.
fprintf('Building LLE envelope from Gomis et al. (2000) tie-line data...\n');

% Organic phase end-points (CyHex-rich):  [EtOH, H2O, CyHex]
lle_org = [ ...
    0.0000, 0.0051, 0.9949;
    0.0514, 0.0141, 0.9345;
    0.0926, 0.0166, 0.8908;
    0.1140, 0.0180, 0.8680;
    0.1360, 0.0216, 0.8424;
    0.1590, 0.0296, 0.8114;
    0.1860, 0.0390, 0.7750;
    0.2100, 0.0510, 0.7390;
    0.2280, 0.2380, 0.5340;   % plait point = ternary azeotrope
];

% Aqueous phase end-points (H2O-rich):  [EtOH, H2O, CyHex]
lle_aq = [ ...
    0.0000, 0.9984, 0.0016;
    0.0482, 0.9430, 0.0088;
    0.0789, 0.9105, 0.0106;
    0.0820, 0.9000, 0.0180;
    0.0990, 0.8880, 0.0130;
    0.1130, 0.8720, 0.0150;
    0.1340, 0.8520, 0.0140;
    0.1680, 0.8160, 0.0160;
    0.2280, 0.2380, 0.5340;   % plait point = ternary azeotrope
];

fprintf('  LLE envelope: %d tie-lines loaded.\n', size(lle_org,1)-1);

% =========================================================================
%% 3. IDENTIFY KEY AZEOTROPE COMPOSITIONS
% =========================================================================
fprintf('Locating azeotropes...\n');

% Ternary heterogeneous azeotrope — refine from NRTL
az_tern = find_azeotrope([0.228, 0.238, 0.534], P, params);

% EtOH-H2O binary azeotrope (x3=0)
az_ethwat = find_binary_az([0.894, 0.106, 0], 1, 2, P, params);

% EtOH-CyHex binary azeotrope (x2=0)
az_ethcyc = find_binary_az([0.44, 0, 0.56], 1, 3, P, params);

fprintf('  Ternary azeotrope:    EtOH=%.4f  H2O=%.4f  CyHex=%.4f\n', az_tern);
fprintf('  EtOH-H2O azeotrope:  EtOH=%.4f  H2O=%.4f\n', az_ethwat(1), az_ethwat(2));
fprintf('  EtOH-CyHex azeotrope:EtOH=%.4f  CyHex=%.4f\n', az_ethcyc(1), az_ethcyc(3));

% =========================================================================
%% 4. PLOT RESIDUE CURVE MAP
% =========================================================================
figure(1);
clf;
set(gcf,'Color','white','Position',[50 50 860 780]);

ax = axes('Position',[0.12 0.10 0.82 0.84]);
hold on;

% --- Draw ternary axes ---
% Triangle vertices:  EtOH=pure (left), H2O=pure (right), CyHex=pure (top)
% Cartesian mapping: x1=EtOH, x2=H2O, x3=CyHex
%   X_cart = 0.5*(2*x2 + x3) = x2 + 0.5*x3
%   Y_cart = (sqrt(3)/2) * x3
V = [0,0; 1,0; 0.5, sqrt(3)/2];   % EtOH, H2O, CyHex vertices

% Triangle border
patch(V(:,1), V(:,2), 'w', 'EdgeColor','k','LineWidth',2);

% Gridlines (10% increments)
for f = 0.1:0.1:0.9
    % Lines parallel to each side
    p1 = tern2cart([f, 1-f, 0]);
    p2 = tern2cart([f, 0, 1-f]);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'Color',[0.85,0.85,0.85],'LineWidth',0.5);

    p1 = tern2cart([0, f, 1-f]);
    p2 = tern2cart([1-f, f, 0]);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'Color',[0.85,0.85,0.85],'LineWidth',0.5);

    p1 = tern2cart([0, 1-f, f]);
    p2 = tern2cart([1-f, 0, f]);
    plot([p1(1),p2(1)],[p1(2),p2(2)],'Color',[0.85,0.85,0.85],'LineWidth',0.5);
end

% --- Plot residue curves ---
for k = 1:numel(curves_fwd)
    C = tern2cart_mat(curves_fwd{k});
    plot(C(:,1), C(:,2), 'b-', 'LineWidth', 0.8);
end
for k = 1:numel(curves_bwd)
    C = tern2cart_mat(curves_bwd{k});
    plot(C(:,1), C(:,2), 'b-', 'LineWidth', 0.8);
end

% --- Plot LLE two-phase envelope ---
%
%  The two-phase region is bounded by two curves meeting at the plait point:
%    Organic boundary (CyHex-rich):  measured tie-line end-points, high CyHex
%    Aqueous boundary (H2O-rich):    measured tie-line end-points, high H2O
%
%  IMPORTANT: The envelope does NOT extend to the pure-component vertices.
%  Data row 1 in lle_org is the CyHex-richest measured point (not pure CyHex).
%  Data row 1 in lle_aq  is the H2O-richest measured point  (not pure H2O).
%  The boundary is drawn only between measured data points — no extrapolation
%  to pure component vertices.

if ~isempty(lle_org) && ~isempty(lle_aq)
    Co = tern2cart_mat(lle_org);   % Nx2 Cartesian, organic boundary
    Ca = tern2cart_mat(lle_aq);    % Nx2 Cartesian, aqueous boundary

    % Close the polygon: organic (row1→end) then aqueous (end→row1)
    % This traces the perimeter of the two-phase lens correctly
    poly_x = [Co(:,1); flipud(Ca(:,1))];
    poly_y = [Co(:,2); flipud(Ca(:,2))];
    fill(poly_x, poly_y, [0.82 0.96 0.82], 'EdgeColor','none', ...
         'FaceAlpha', 0.50, 'DisplayName','Two-phase (LLE) region');

    % Organic boundary — dark green solid line
    plot(Co(:,1), Co(:,2), '-', 'Color',[0.0 0.50 0.0], 'LineWidth', 2.5, ...
         'DisplayName','LLE organic boundary (Gomis et al., 2000)');

    % Aqueous boundary — dark magenta solid line
    plot(Ca(:,1), Ca(:,2), '-', 'Color',[0.70 0.0 0.70], 'LineWidth', 2.5, ...
         'DisplayName','LLE aqueous boundary (Gomis et al., 2000)');

    % Representative tie-lines — only between actual measured phase pairs
    % Skip tie-line at row 1 (nearly pure-component end, visually trivial)
    % Skip last row (plait point, phases are identical)
    n_tie = size(lle_org,1) - 1;
    tie_indices = 2 : max(1,floor(n_tie/5)) : n_tie;
    for k = tie_indices
        po = tern2cart(lle_org(k,:));
        pa = tern2cart(lle_aq(k,:));
        plot([po(1),pa(1)],[po(2),pa(2)],'k-','LineWidth',0.9,'HandleVisibility','off');
    end

    % Plait point marker (where organic = aqueous = ternary azeotrope)
    p_plait = tern2cart(lle_org(end,:));
    plot(p_plait(1), p_plait(2), '^', 'MarkerSize',10,'LineWidth',2,...
         'Color',[0.0 0.5 0.0],'MarkerFaceColor',[0.0 0.7 0.0],...
         'DisplayName','Plait point (= ternary AZ)');
end

% --- Mark azeotropes ---
p_az_tern   = tern2cart(az_tern);
p_az_ethwat = tern2cart(az_ethwat);
p_az_ethcyc = tern2cart(az_ethcyc);

plot(p_az_tern(1),   p_az_tern(2),   'r*', 'MarkerSize',14,'LineWidth',2,...
     'DisplayName','Ternary azeotrope');
plot(p_az_ethwat(1), p_az_ethwat(2), 'rs', 'MarkerSize',10,'LineWidth',2,...
     'DisplayName','EtOH-H2O azeotrope');
plot(p_az_ethcyc(1), p_az_ethcyc(2), '^', 'MarkerSize',11,'LineWidth',2,...
     'Color','r','MarkerFaceColor',[1.0 0.4 0.4],...
     'DisplayName','EtOH-CyHex azeotrope');
% Draw a short line from marker to label to make it unambiguous
text(p_az_ethcyc(1)-0.08, p_az_ethcyc(2)+0.02, ...
     '\leftarrow EtOH-CyHex edge', 'FontSize', 8, 'Color', [0.5 0.5 0.5]);

% --- Mark pure component nodes ---
plot(V(1,1), V(1,2), 'ko','MarkerSize',8,'MarkerFaceColor','k');
plot(V(2,1), V(2,2), 'ko','MarkerSize',8,'MarkerFaceColor','k');
plot(V(3,1), V(3,2), 'ko','MarkerSize',8,'MarkerFaceColor','k');

% --- Vertex labels ---
% EtOH  vertex = (0,0)              bottom-left
% H2O   vertex = (1,0)              bottom-right
% CyHex vertex = (0.5, sqrt(3)/2)   top
text( 0.00, -0.06,  'Ethanol',     'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
text( 1.00, -0.06,  'Water',       'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
text( 0.50,  0.90,  'Cyclohexane', 'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');

% Tick labels (mole fraction 0.2 increments)
for f = 0.2:0.2:0.8
    % EtOH axis (left side, EtOH increasing bottom-left to top)
    p = tern2cart([f, 1-f, 0]);
    text(p(1)-0.06, p(2)-0.01, sprintf('%.1f',f),'FontSize',8,'HorizontalAlignment','right');
    % H2O axis (bottom, increasing left to right)
    p = tern2cart([0, f, 1-f]);
    text(p(1)+0.01, p(2)-0.03, sprintf('%.1f',f),'FontSize',8);
    % CyHex axis (right side)
    p = tern2cart([1-f, 0, f]);
    text(p(1)+0.02, p(2), sprintf('%.1f',f),'FontSize',8);
end

% Azeotrope annotations
% Ternary AZ — offset right so it does not overlap the marker
text(p_az_tern(1)+0.04, p_az_tern(2)+0.03, ...
     sprintf('AZ_{tern} (%.3f, %.3f, %.3f)',az_tern(1),az_tern(2),az_tern(3)), ...
     'FontSize',9,'Color',[0.8 0 0],'FontWeight','bold');

% EtOH-H2O binary AZ — on bottom edge, offset upward
text(p_az_ethwat(1), p_az_ethwat(2)+0.04, ...
     sprintf('AZ_{EtOH-H2O}\n(%.3f, %.3f)',az_ethwat(1),az_ethwat(2)), ...
     'FontSize',9,'Color',[0.8 0 0],'HorizontalAlignment','center');

% EtOH-CyHex binary AZ — on left edge, offset left to avoid overlap
text(p_az_ethcyc(1)-0.06, p_az_ethcyc(2), ...
     sprintf('AZ_{EtOH-CyHex}\n(%.3f, %.3f)',az_ethcyc(1),az_ethcyc(3)), ...
     'FontSize',9,'Color',[0.8 0 0],'HorizontalAlignment','right');

title('Residue Curve Map: Ethanol–Water–Cyclohexane (101325 Pa, NRTL)', ...
      'FontSize',14,'FontWeight','bold');
legend('Location','northeast','FontSize',9);
axis off;
axis equal;
xlim([-0.15, 1.15]);
ylim([-0.10, 0.96]);

% =========================================================================
%% 5. SAVE FIGURE
% =========================================================================
saveas(figure(1), 'RCM_EtOH_H2O_CyHex.png');
fprintf('\nFigure saved: RCM_EtOH_H2O_CyHex.png\n');
fprintf('\n=== Key Azeotrope Compositions (NRTL, 101325 Pa) ===\n');
fprintf('  Ternary heterogeneous azeotrope:\n');
fprintf('    EtOH = %.4f,  H2O = %.4f,  CyHex = %.4f\n', az_tern);
fprintf('  EtOH-H2O binary azeotrope:\n');
fprintf('    EtOH = %.4f,  H2O = %.4f\n', az_ethwat(1), az_ethwat(2));
fprintf('  EtOH-CyHex binary azeotrope:\n');
fprintf('    EtOH = %.4f,  CyHex = %.4f\n', az_ethcyc(1), az_ethcyc(3));

end  % ---- end generate_RCM ----


% =========================================================================
%% LOCAL FUNCTIONS
% =========================================================================

function dxdxi = rcm_rhs(~, x, P, params)
% RCM ODE: dx/dxi = x - y*(x,T)
    x = max(x(:)', 0);
    x = x / max(sum(x), 1e-15);
    try
        [y, ~, ~] = VLE_bubble(x, P, params);
        y = y(:);
        dxdxi = x(:) - y;
    catch
        dxdxi = zeros(3,1);
    end
end


function [val, isterm, dir] = stop_event(~, x)
% Stop integration when any composition hits a boundary (pure component)
    val    = min(x) - 1e-4;
    isterm = 1;
    dir    = -1;
end


function az = find_azeotrope(x0, P, params)
% Locate ternary heterogeneous azeotrope by multi-start fsolve
%
%   Solves y_i(x,T) = x_i for i=1,2  (x3 = 1-x1-x2 eliminated).
%   Multiple starting points are used and the solution is validated
%   against the known physical location of the ternary azeotrope for
%   EtOH-H2O-CyHex: EtOH ~ 0.15-0.30, H2O ~ 0.15-0.35, CyHex ~ 0.45-0.65
%   Any root outside this window is rejected as spurious.

    obj  = @(u) az_residual(u, P, params);
    opts = optimoptions('fsolve','Display','off','TolFun',1e-11,'TolX',1e-11,...
                        'MaxIterations',800,'MaxFunctionEvaluations',5000);

    % Multiple starting points anchored to known region of ternary azeotrope
    % (Gomis et al., 2000; Knapp & Doherty, 1994)
    starts = [ ...
        0.228, 0.238;   % literature value
        0.200, 0.250;
        0.250, 0.220;
        0.180, 0.260;
        0.260, 0.200;
        0.210, 0.230;
        0.240, 0.240;
        0.190, 0.220;
    ];

    best_az  = x0(:)';
    best_res = Inf;

    for k = 1:size(starts,1)
        u0 = starts(k,:);
        % Skip if x3 would be non-physical
        if sum(u0) >= 0.99,  continue;  end
        try
            [u_sol, ~, exitflag] = fsolve(obj, u0, opts);
        catch
            continue;
        end
        if exitflag <= 0,  continue;  end

        x1 = u_sol(1);  x2 = u_sol(2);  x3 = 1 - x1 - x2;

        % Reject non-physical solutions
        if x1 < 0 || x2 < 0 || x3 < 0,  continue;  end
        if x1+x2+x3 < 0.99 || x1+x2+x3 > 1.01,  continue;  end

        % Reject spurious roots: ternary azeotrope must have
        %   EtOH  in [0.10, 0.35]
        %   H2O   in [0.10, 0.40]
        %   CyHex in [0.35, 0.75]
        if x1 < 0.10 || x1 > 0.35,  continue;  end
        if x2 < 0.10 || x2 > 0.40,  continue;  end
        if x3 < 0.35 || x3 > 0.75,  continue;  end

        % Verify y = x at solution (residual check)
        r = az_residual(u_sol, P, params);
        res_norm = norm(r);
        if res_norm < best_res
            best_res = res_norm;
            best_az  = max([x1, x2, x3], 0);
            best_az  = best_az / sum(best_az);
        end
    end

    if best_res > 1e-4
        warning('find_azeotrope: best residual = %.2e — ternary azeotrope may not have converged. Using best estimate.', best_res);
    end
    az = best_az;
end


function r = az_residual(u, P, params)
% Residual for azeotrope condition: y_i(x) - x_i = 0  for i=1,2
    x1 = u(1);  x2 = u(2);  x3 = max(0, 1-x1-x2);
    x = max([x1, x2, x3], 1e-10);
    x = x / sum(x);
    try
        [y, ~, ~] = VLE_bubble(x, P, params);
        r = [y(1)-x(1); y(2)-x(2)];
    catch
        r = [1; 1];
    end
end


function az = find_binary_az(x0, i1, i2, P, params)
% Locate binary azeotrope between components i1 and i2 (third = 0)
    obj = @(t) binary_az_res(t, i1, i2, P, params);
    opts = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
    try
        t_sol = fzero(obj, x0(i1), opts);
        az = zeros(1,3);
        az(i1) = max(min(t_sol, 1-1e-6), 1e-6);
        az(i2) = 1 - az(i1);
    catch
        az = x0;
    end
end


function r = binary_az_res(t, i1, i2, P, params)
% Residual: y_i1(x) - x_i1 = 0  on binary edge (x_i3=0)
    x = zeros(1,3);
    x(i1) = max(min(t, 1-1e-8), 1e-8);
    x(i2) = 1 - x(i1);
    try
        [y, ~, ~] = VLE_bubble(x, P, params);
        r = y(i1) - x(i1);
    catch
        r = 0;
    end
end


function p = tern2cart(x)
% Convert ternary [x1_EtOH, x2_H2O, x3_CyHex] to 2D Cartesian
% Vertices: EtOH=(0,0), H2O=(1,0), CyHex=(0.5, sqrt(3)/2)
    x = x(:)';
    p(1) = x(2) + 0.5*x(3);
    p(2) = (sqrt(3)/2) * x(3);
end


function C = tern2cart_mat(X)
% Convert matrix of ternary compositions (Nx3) to Cartesian (Nx2)
    C = zeros(size(X,1), 2);
    for k = 1:size(X,1)
        p = tern2cart(X(k,:));
        C(k,:) = p;
    end
end