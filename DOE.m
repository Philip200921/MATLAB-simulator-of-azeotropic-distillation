function validate_DOE_solutions()
clc;

%  FIXED BASE PARAMETERS

F        = 460;
conv     = 0.9999;
P1       = 101325;   P2       = 101325;
n2       = 25;       R2       = 0.30;
ethfracb = 0.9999;   H2Ofracb = 0.9999;   ethfracd = 0.10;

DE_solutions = [ ...
    1,  3.223,  0.100,  76.008,  8.897,  11200.748,  0.870;
    2,  3.223,  0.100,  75.692,  8.897,  11200.822,  0.870;
    3,  3.223,  0.100,  76.454,  8.897,  11200.995,  0.870;
    4,  3.223,  0.100,  76.810,  8.897,  11201.930,  0.870;
    5,  3.224,  0.100,  74.840,  8.897,  11203.494,  0.870  ];

N_sol = size(DE_solutions, 1);

%  PREALLOCATE SIMULATION RESULTS

n1_used    = nan(N_sol, 1);   % actual integer stages used
D1_all     = nan(N_sol, 1);   % distillate flow     [kmol/hr]
MFED_sim   = nan(N_sol, 1);   % simulated MFED      [kmol/hr]
Energy_sim = nan(N_sol, 1);   % simulated energy    [kW]
err_MFED   = nan(N_sol, 1);   % absolute error MFED [kmol/hr]
err_Energy = nan(N_sol, 1);   % absolute error duty [kW]
pct_MFED   = nan(N_sol, 1);   % relative error MFED [%]
pct_Energy = nan(N_sol, 1);   % relative error duty [%]
status     = cell(N_sol, 1);

%  SIMULATION LOOP

for k = 1:N_sol
    sol_no  = DE_solutions(k, 1);
    R1_k    = DE_solutions(k, 2);
    SFR_k   = DE_solutions(k, 3);
    n1_raw  = DE_solutions(k, 4);
    n1_k    = round(n1_raw);      % runsimulation requires integer stages

    % SFR to feed mole fractions
    ent_k = SFR_k / (1 + SFR_k);
    eth_k = 0.8333 * (1 - ent_k);
    wat_k = 0.1667 * (1 - ent_k);

    % D1 from material balance
    B1_k = (conv * F * eth_k) / ethfracb;
    D1_k = F - B1_k;

    n1_used(k) = n1_k;
    D1_all(k)  = D1_k;

   % fprintf('  Solution %d | R1=%.3f  SFR=%.3f  Stages(raw)=%.3f → n1=%d\n', ...
   %          sol_no, R1_k, SFR_k, n1_raw, n1_k);
   %  fprintf('             | ent=%.5f  eth=%.5f  wat=%.5f  D1=%.4f kmol/hr\n', ...
   %          ent_k, eth_k, wat_k, D1_k);

    try
        warning('off', 'all');
        [Res, ~] = runsimulation(F, conv, P1, P2, n1_k, n2, ...
                                  R1_k, R2, ...
                                  eth_k, wat_k, ent_k, ...
                                  ethfracb, H2Ofracb, ethfracd);
        warning('on', 'all');

        xTop_eth = Res.Column1.xTop(1);
        Qreb     = Res.Column1.Qreb_kW;
        Qcond    = Res.Column1.Qcond_kW;

        if ~isfinite(xTop_eth) || ~isfinite(Qreb) || ~isfinite(Qcond)
            error('Non-finite output from runsimulation.');
        end

        MFED_sim(k)   = D1_k * xTop_eth;
        Energy_sim(k) = abs(Qreb) + abs(Qcond);
        err_MFED(k)   = MFED_sim(k)   - DE_solutions(k, 5);
        err_Energy(k) = Energy_sim(k)  - DE_solutions(k, 6);
        pct_MFED(k)   = abs(err_MFED(k))   / DE_solutions(k, 5) * 100;
        pct_Energy(k) = abs(err_Energy(k))  / DE_solutions(k, 6) * 100;
        status{k}     = 'OK';
        % 
        % fprintf('             | MFED_sim=%.4f kmol/hr  Energy_sim=%.2f kW\n\n', ...
        %         MFED_sim(k), Energy_sim(k));

    catch ME
        status{k} = ME.message;
        fprintf('             | FAILED — %s\n\n', ME.message);
    end
end

% Comparison talble

fprintf('  (Simulation vs Design-Expert Response Surface)\n');


fprintf('%-4s %-5s %-5s %-6s %-4s | %-10s %-10s %-8s %-6s | %-12s %-12s %-8s %-6s\n', ...
        'No.','R1','SFR','Stages','n1', ...
        'MFED_DE','MFED_Sim','Err','Err%', ...
        'Energy_DE','Energy_Sim','Err','Err%');
fprintf('%s\n', repmat('-',1,110));

for k = 1:N_sol
    if strcmp(status{k}, 'OK')
        fprintf('%-4d %-5.3f %-5.3f %-6.3f %-4d | %-10.4f %-10.4f %-8.4f %-6.2f | %-12.2f %-12.2f %-8.2f %-6.2f\n', ...
                DE_solutions(k,1), DE_solutions(k,2), DE_solutions(k,3), ...
                DE_solutions(k,4), n1_used(k), ...
                DE_solutions(k,5), MFED_sim(k), err_MFED(k), pct_MFED(k), ...
                DE_solutions(k,6), Energy_sim(k), err_Energy(k), pct_Energy(k));
    else
        fprintf('%-4d %-5.3f %-5.3f %-6.3f %-4d | %-10.4f %-10s %-8s %-6s | %-12.2f %-12s %-8s %-6s  FAILED\n', ...
                DE_solutions(k,1), DE_solutions(k,2), DE_solutions(k,3), ...
                DE_solutions(k,4), n1_used(k), ...
                DE_solutions(k,5), '---', '---', '---', ...
                DE_solutions(k,6), '---', '---', '---');
    end
end
fprintf('%s\n', repmat('-',1,110));

% Summary statistics
ok_idx = strcmp(status, 'OK');
if any(ok_idx)
    fprintf('\nMean absolute error — MFED  : %.4f kmol/hr  (%.2f%%)\n', ...
            mean(abs(err_MFED(ok_idx))), mean(pct_MFED(ok_idx)));
    fprintf('Mean absolute error — Energy: %.2f kW       (%.2f%%)\n', ...
            mean(abs(err_Energy(ok_idx))), mean(pct_Energy(ok_idx)));
end