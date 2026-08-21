function run_DOE_table()
%   Factor levels:
%     A — Reflux ratio R1 : 0.5 (low),  2.0 (mid),  3.5 (high)
%     B — SFR             : 0.11 (low), 1.55 (mid), 3.00 (high)
%     C — Stages n1       : 10 (low),   45 (mid),   80 (high)
%
%   SFR-to-feed-composition conversion (EtOH:H2O ratio preserved at 5:1):
%     ent = SFR / (1 + SFR)
%     eth = 0.8333 * (1 - ent)
%     wat = 0.1667 * (1 - ent)

clc;
F        = 460;
conv     = 0.9999;
P1       = 101325;   P2       = 101325;
n2       = 25;       R2       = 0.30;
ethfracb = 0.9999;   H2Ofracb = 0.9999;   ethfracd = 0.10;

%         Std  Run   R1    SFR    n1
DOE = [ ...
         6,   1,   3.5,  1.55,  10;
        11,   2,   2.0,  0.11,  80;
         1,   3,   0.5,  0.11,  45;
         8,   4,   3.5,  1.55,  80;
        12,   5,   2.0,  3.00,  80;
         2,   6,   3.5,  0.11,  45;
        10,   7,   2.0,  3.00,  10;
         7,   8,   0.5,  1.55,  80;
        13,   9,   2.0,  1.55,  45;
         5,  10,   0.5,  1.55,  10;
         4,  11,   3.5,  3.00,  45;
         9,  12,   2.0,  0.11,  10;
         3,  13,   0.5,  3.00,  45  ];

N_runs = size(DOE, 1);

%  PREALLOCATE

MFED_pct = nan(N_runs, 1);
Q_total  = nan(N_runs, 1);
status   = cell(N_runs, 1);

%  SIMULATION LOOP

fprintf('Running %d DOE simulations...\n\n', N_runs);
fprintf('%-4s %-4s %-5s %-6s %-7s | %-8s %-8s %-8s | %-12s %-12s\n', ...
        'Std','Run','R1','SFR','Stages','ent','eth','wat','MFED (%)','Energy (kW)');
fprintf('%s\n', repmat('-', 1, 85));

for k = 1:N_runs
    std_k = DOE(k, 1);
    run_k = DOE(k, 2);
    R1_k  = DOE(k, 3);
    SFR_k = DOE(k, 4);
    n1_k  = DOE(k, 5);

    % SFR to feed mole fractions
    ent_k = SFR_k / (1 + SFR_k);
    eth_k = 0.8333 * (1 - ent_k);
    wat_k = 0.1667 * (1 - ent_k);
    % 
    fprintf('%-4d %-4d %-5.1f %-6.2f %-7d | %-8.4f %-8.4f %-8.4f | ', ...
            std_k, run_k, R1_k, SFR_k, n1_k, ent_k, eth_k, wat_k);

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
            error('Non-finite output.');
        end

        MFED_pct(k) = xTop_eth * 100;
        Q_total(k)  = abs(Qreb) + abs(Qcond);
        status{k}   = 'OK';

        fprintf('%-12.4f %-12.2f  OK\n', MFED_pct(k), Q_total(k));

    catch ME
        status{k} = ME.message;
        fprintf('%-12s %-12s  FAILED — %s\n', '---', '---', ME.message);
    end
end

%  FINAL TABLE

fprintf('DOE RESPONSE TABLE\n');


fprintf('%-4s  %-4s  %-12s  %-8s  %-8s  |  %-12s  %-12s\n', ...
        'Std', 'Run', 'A:Refl.Ratio', 'B:SFR', 'C:Stages', ...
        'MFED (%)', 'Energy (kW)');
fprintf('%s\n', repmat('-', 1, 72));

for k = 1:N_runs
    if strcmp(status{k}, 'OK')
        fprintf('%-4d  %-4d  %-12.1f  %-8.2f  %-8d  |  %-12.4f  %-12.2f\n', ...
                DOE(k,1), DOE(k,2), DOE(k,3), DOE(k,4), DOE(k,5), ...
                MFED_pct(k), Q_total(k));
    else
        fprintf('%-4d  %-4d  %-12.1f  %-8.2f  %-8d  |  %-12s  %-12s  [FAILED]\n', ...
                DOE(k,1), DOE(k,2), DOE(k,3), DOE(k,4), DOE(k,5), ...
                '---', '---');
    end
end
fprintf('%s\n', repmat('-', 1, 72));

n_ok     = sum(strcmp(status, 'OK'));
n_failed = N_runs - n_ok;
fprintf('\nSuccessful: %d / %d    Failed: %d / %d\n', n_ok, N_runs, n_failed, N_runs);
end