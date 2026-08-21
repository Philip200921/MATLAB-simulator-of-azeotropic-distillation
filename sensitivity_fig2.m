function sensitivity_fig2()
% sensitivity_fig_total_duty
%   Effect of number of stages (n1) and reflux ratio (R1) on
%   TOTAL Column 1 heat duty:
%
%       Q_total = Q_reboiler + |Q_condenser|   [kW]
%
%   Six R1 values are plotted as separate curves with n1 on the x-axis.
%   Each (R1, n1) point calls runsimulation exactly ONCE, and both
%   Qreb and Qcond are extracted from the single result structure.
%
%   USAGE:  >> sensitivity_fig_total_duty

clc;
fprintf('=============================================\n');
fprintf('  Total Heat Duty: n1 and R1 sensitivity\n');
fprintf('=============================================\n\n');

% =========================================================================
%  BASE-CASE PARAMETERS
% =========================================================================
F        = 170;
conv     = 0.9999;
P1       = 101325;   P2       = 101325;
n2       = 25;       R2       = 0.30;
ethfracf = 0.50;     watfracf = 0.10;   entfracf = 0.40;
ethfracb = 0.9999;   H2Ofracb = 0.9999; ethfracd = 0.10;

% =========================================================================
%  SWEEP VECTORS
% =========================================================================
n1_values = round(linspace(10, 80, 12));
R1_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];

n_R1 = numel(R1_values);
n_n1 = numel(n1_values);

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
%  SIMULATION LOOP  (single call per point — both Qreb and Qcond extracted)
% =========================================================================
Q_total = nan(n_R1, n_n1);   % Q_reboiler + |Q_condenser|  [kW]

for r = 1:n_R1
    for k = 1:n_n1
        fprintf('  R1 = %.1f  |  n1 = %2d  ...  ', R1_values(r), n1_values(k));
        try
            [Res, ~] = runsimulation(F, conv, P1, P2, ...
                                     n1_values(k), n2, ...
                                     R1_values(r), R2, ...
                                     ethfracf, watfracf, entfracf, ...
                                     ethfracb, H2Ofracb, ethfracd);

            Qreb  = Res.Column1.Qreb_kW;
            Qcond = Res.Column1.Qcond_kW;

            % Validate both outputs before accepting the point
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
fig = figure(1); clf;
set(fig, 'Color', 'white', 'Position', [200 100 720 520]);

hold on;
for r = 1:n_R1
    % Only plot points where the simulation succeeded
    valid = isfinite(Q_total(r, :));
    if ~any(valid)
        warning('sensitivity_fig_total_duty: no valid points for R1 = %.1f — curve skipped.', ...
                R1_values(r));
        continue;
    end
    plot(n1_values(valid), Q_total(r, valid), ...
         ['-' markers{r}], ...
         'Color',           colors(r, :), ...
         'LineWidth',        1.8, ...
         'MarkerSize',       7, ...
         'MarkerFaceColor',  colors(r, :), ...
         'DisplayName',      sprintf('R_1 = %.1f', R1_values(r)));
end
hold off;

xlabel('Number of theoretical stages,  n_1  (–)', ...
       'FontSize', 12, 'FontWeight', 'bold');
ylabel('Total heat duty,  Q_{total} = Q_{reb} + |Q_{cond}|  (kW)', ...
       'FontSize', 11, 'FontWeight', 'bold');
title('Effect of number of stages and reflux ratio on total Column 1 heat duty', ...
      'FontSize', 11, 'FontWeight', 'bold');

legend('Location', 'northeast', 'FontSize', 9, 'Box', 'on');
grid on;
box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0);

% Export publication-quality files
exportgraphics(fig, 'Fig_TotalDuty_n1_R1.pdf', 'ContentType', 'vector', 'Resolution', 300);
exportgraphics(fig, 'Fig_TotalDuty_n1_R1.png', 'Resolution', 300);
fprintf('\nFigure saved: Fig_TotalDuty_n1_R1.pdf / .png\n');
end