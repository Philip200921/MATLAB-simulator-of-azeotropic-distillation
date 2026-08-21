function sensitivity_fig3()

%   Six R1 values are plotted as separate curves with S/F on the x-axis.
%   Entrainer mole fraction swept from 0.10 to 0.75 (S/F: 0.111 to 3.000).
%   Each (R1, S/F) point calls runsimulation exactly ONCE.


clc;
fprintf('=============================================\n');
fprintf('  Total Heat Duty: S/F and R1 sensitivity\n');
fprintf('=============================================\n\n');

% =========================================================================
%  BASE-CASE PARAMETERS
% =========================================================================
F        = 460;
conv     = 0.9999;
P1       = 101325;   P2       = 101325;
n1       = 30;       n2       = 25;       R2       = 0.30;
eth0     = 0.50;     wat0     = 0.10;
ethfracb = 0.9999;   H2Ofracb = 0.9999;   ethfracd = 0.10;

% =========================================================================
%  SWEEP VECTORS
% =========================================================================
ent_values = linspace(0.10, 0.75, 15);   % entrainer mole fraction sweep
R1_values  = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];

n_R1  = numel(R1_values);
n_ent = numel(ent_values);

% =========================================================================
%  PLOT STYLE
% =========================================================================
colors  = [0.85 0.10 0.10;
           1.00 0.55 0.00;
           0.10 0.60 0.10;
           0.10 0.40 0.85;
           0.50 0.00 0.50;
           0.00 0.00 0.00];

markers = {'o', 's', 'd', '^', 'v', 'p'};

% =========================================================================
%  SIMULATION LOOP  (single call per point)
% =========================================================================
Q_total = nan(n_R1, n_ent);

for r = 1:n_R1
    for k = 1:n_ent
        ent_k = ent_values(k);

        % Scale EtOH and H2O proportionally so all fractions sum to 1
        eth_k = eth0 * (1 - ent_k) / (eth0 + wat0);
        wat_k = wat0 * (1 - ent_k) / (eth0 + wat0);

        fprintf('  R1 = %.1f  |  S/F = %.3f  (EtOH=%.3f  H2O=%.3f  CyHex=%.3f)  ...  ', ...
                R1_values(r), ent_k/(1-ent_k), eth_k, wat_k, ent_k);
        try
            [Res, ~] = runsimulation(F, conv, P1, P2, n1, n2, ...
                                     R1_values(r), R2, ...
                                     eth_k, wat_k, ent_k, ...
                                     ethfracb, H2Ofracb, ethfracd);

            Qreb  = Res.Column1.Qreb_kW;
            Qcond = Res.Column1.Qcond_kW;

            if isfinite(Qreb) && isfinite(Qcond)
                Q_total(r, k) = abs(Qreb) + abs(Qcond);
                fprintf('Q_total = %.2f kW\n', Q_total(r, k));
            else
                fprintf('non-finite duty — skipped\n');
            end

        catch ME
            fprintf('FAILED — %s\n', ME.message);
        end
    end
end

% =========================================================================
%  PLOT
% =========================================================================
SF_ratio = ent_values ./ (1 - ent_values);   % convert entrainer fraction to S/F

fig = figure(1); clf;
set(fig, 'Color', 'white', 'Position', [250 100 720 520]);

hold on;
for r = 1:n_R1
    valid = isfinite(Q_total(r, :));
    if ~any(valid)
        warning('sensitivity_fig3: no valid points for R1 = %.1f — curve skipped.', ...
                R1_values(r));
        continue;
    end
    plot(SF_ratio(valid), Q_total(r, valid), ...
         ['-' markers{r}], ...
         'Color',           colors(r, :), ...
         'LineWidth',        1.8, ...
         'MarkerSize',       7, ...
         'MarkerFaceColor',  colors(r, :), ...
         'DisplayName',      sprintf('R_1 = %.1f', R1_values(r)));
end
hold off;

xlabel('Solvent-to-feed ratio,  S/F  (–)', ...
       'FontSize', 12, 'FontWeight', 'bold');
ylabel('Total heat duty,  Q_{total} = Q_{reb} + |Q_{cond}|  (kW)', ...
       'FontSize', 11, 'FontWeight', 'bold');
title('Effect of solvent flowrate and reflux ratio on total Column 1 heat duty', ...
      'FontSize', 11, 'FontWeight', 'bold');

legend('Location', 'northwest', 'FontSize', 9, 'Box', 'on');
grid on;
box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0);
end