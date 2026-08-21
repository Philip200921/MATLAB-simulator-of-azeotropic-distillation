function sensitivity_fig1()
% sensitivity_fig1  Figure 1 — Impact of R1 on distillate composition and heat duty
%
%   Dual y-axis: left = EtOH mole fraction in Col1 top vapour
%                right = Reboiler and Condenser duty (kW)
%   Sweep R1 from 0.3 to 3.5 in 12 steps.
%
%   USAGE:  >> sensitivity_fig1

clear functions
clc;
fprintf('=============================================\n');
fprintf('  Figure 1: R1 vs Distillate Composition\n');
fprintf('           and Heat Duty\n');
fprintf('=============================================\n\n');

% Base-case parameters
F        = 170;   conv     = 0.9999;
P1       = 101325; P2      = 101325;
n1       = 30;    n2       = 25;
R2       = 0.30;
ethfracf = 0.50;  watfracf = 0.10;  entfracf = 0.40;
ethfracb = 0.9999; H2Ofracb = 0.9999; ethfracd = 0.10;

R1_values = linspace(0.3, 3.5, 12);
xD_eth  = nan(size(R1_values));
xD_wat  = nan(size(R1_values));
xD_ent  = nan(size(R1_values));
Qreb1   = nan(size(R1_values));
Qcond1  = nan(size(R1_values));

for k = 1:numel(R1_values)
    fprintf('  R1 = %.3f ... ', R1_values(k));
    try
        [Res, ~] = runsimulation(F, conv, P1, P2, n1, n2, R1_values(k), R2, ...
                                 ethfracf, watfracf, entfracf, ...
                                 ethfracb, H2Ofracb, ethfracd);
        xD_eth(k)  = Res.Column1.xTop(1);
        xD_wat(k)  = Res.Column1.xTop(2);
        xD_ent(k)  = Res.Column1.xTop(3);
        Qreb1(k)   = abs(Res.Column1.Qreb_kW);
        Qcond1(k)  = abs(Res.Column1.Qcond_kW);
        fprintf('xD_EtOH = %.4f | Qreb = %.1f kW\n', xD_eth(k), Qreb1(k));
    catch ME
        fprintf('FAILED — %s\n', ME.message);
    end
end

% Plot — dual y-axis
fig = figure(1); clf;
set(fig, 'Color', 'white', 'Position', [100 100 680 500]);

yyaxis left
plot(R1_values, xD_eth, '-o', 'Color', [0.85 0.10 0.10], ...
     'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.10 0.10], ...
     'DisplayName', 'EtOH in distillate');
hold on;
plot(R1_values, xD_wat, '-s', 'Color', [0.10 0.40 0.85], ...
     'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', [0.10 0.40 0.85], ...
     'DisplayName', 'H_2O in distillate');
plot(R1_values, xD_ent, '-d', 'Color', [0.50 0.00 0.50], ...
     'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', [0.50 0.00 0.50], ...
     'DisplayName', 'CyHex in distillate');
ylabel('Mole fraction in distillate', 'FontSize', 11, 'Color', 'k');
set(gca, 'YColor', 'k');

yyaxis right
plot(R1_values, Qreb1,  '--^', 'Color', [0.10 0.60 0.10], ...
     'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', [0.10 0.60 0.10], ...
     'DisplayName', 'Reboiler duty Q_r');
plot(R1_values, Qcond1, '--v', 'Color', [1.00 0.55 0.00], ...
     'LineWidth', 1.8, 'MarkerSize', 7, 'MarkerFaceColor', [1.00 0.55 0.00], ...
     'DisplayName', 'Condenser duty Q_c');
ylabel('Energy duty (kW)', 'FontSize', 11, 'Color', 'k');
set(gca, 'YColor', 'k');

xlabel('Reflux ratio  R_1', 'FontSize', 12, 'FontWeight', 'bold');
title('Figure 1.  Impact of reflux ratio on distillate composition and heat duty', ...
      'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 9);
grid on; box on;
set(gca, 'FontSize', 10, 'LineWidth', 1.0, 'XLim', [0.2 3.6]);
end