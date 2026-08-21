function gamma = activity_coeff_NRTL(x, T, params)
% activity_coeff_NRTL  Activity coefficients via the NRTL model
%
%   gamma = activity_coeff_NRTL(x, T, params)
%
%   Implements the multicomponent NRTL equation (Renon & Prausnitz, 1968):
%
%       ln(gamma_i) = [ sum_j(x_j*tau_ji*G_ji) ] / [ sum_k(x_k*G_ki) ]
%                   + sum_j { [ x_j*G_ij / sum_k(x_k*G_kj) ]
%                             * [ tau_ij - sum_k(x_k*tau_kj*G_kj)
%                                          / sum_k(x_k*G_kj) ] }
%
%   where:
%       tau_ij = Aij + Bij/T           (binary interaction parameter)
%       G_ij   = exp(-alpha_ij * tau_ij)
%       alpha  = non-randomness parameter (symmetric: alpha_ij = alpha_ji)
%       tau_ii = 0,  G_ii = 1  (by definition)
%
%   Reference: Renon H., Prausnitz J.M. (1968) AIChE J., 14(1), 135-144.
%
%   Input:
%       x      – 1xN liquid mole fraction row vector (N = number of components)
%       T      – scalar temperature [K]
%       params – structure with fields:
%                  .NRTL.Aij    NxN matrix of A_ij coefficients
%                  .NRTL.Bij    NxN matrix of B_ij coefficients
%                  .NRTL.alpha  NxN symmetric non-randomness matrix
%
%   Output:
%       gamma  – 1xN row vector of activity coefficients (dimensionless)
%
%   Notes:
%     * Diagonal elements of Aij, Bij must be zero (tau_ii = 0 → G_ii = 1).
%     * For publishable work, validate gamma against DECHEMA VLE data.

    % --- input validation ---------------------------------------------------
    x = x(:)';          % enforce row vector
    N = numel(x);

    if ~isscalar(T) || T <= 0
        error('activity_coeff_NRTL: T must be a positive scalar [K].');
    end
    if abs(sum(x) - 1) > 1e-6
        warning('activity_coeff_NRTL: mole fractions sum to %.6f (not 1); normalising.', sum(x));
        x = x / sum(x);
    end

    x = max(x, 0);   % prevent tiny negative values from propagating

    % --- retrieve parameters ------------------------------------------------
    Aij   = params.NRTL.Aij;     % NxN
    Bij   = params.NRTL.Bij;     % NxN
    alpha = params.NRTL.alpha;   % NxN  (non-randomness; symmetric)

    % --- compute tau and G matrices -----------------------------------------
    Tau = Aij + Bij ./ T;              % NxN  – tau_ij(T)
    G   = exp(-alpha .* Tau);          % NxN  – G_ij(T)

    % Enforce diagonal convention (tau_ii = 0 → G_ii = 1)
    for i = 1:N
        Tau(i,i) = 0;
        G(i,i)   = 1;
    end

    % --- compute ln(gamma_i) for each component i --------------------------
    lngamma = zeros(1, N);

    for i = 1:N
        % Term 1:  sum_j(x_j * tau_ji * G_ji) / sum_k(x_k * G_ki)
        denom1  = max(x * G(:,i), 1e-300);   % sum_k(x_k * G_ki)  [scalar]
        numer1  = x * (Tau(:,i) .* G(:,i));  % sum_j(x_j*tau_ji*G_ji) [scalar]
        term1   = numer1 / denom1;

        % Term 2:  sum_j { x_j*G_ij/[sum_k x_k*G_kj] * [tau_ij - T2_inner_j] }
        term2 = 0;
        for j = 1:N
            denom2  = max(x * G(:,j), 1e-300);            % sum_k(x_k*G_kj)
            numer2  = x * (Tau(:,j) .* G(:,j));           % sum_k(x_k*tau_kj*G_kj)
            t2inner = numer2 / denom2;
            term2   = term2 + x(j) * G(i,j) / denom2 * (Tau(i,j) - t2inner);
        end

        lngamma(i) = term1 + term2;
    end

    gamma = exp(lngamma);
end














% function gamma = activity_coeff_NRTL(x, T, params_local)
%     % Compute activity coefficients using NRTL model
%     R = 8.314;  % J/mol-K
%     Aij = params_local.NRTL.Aij;
%     Bij = params_local.NRTL.Bij;
%     Cij = params_local.NRTL.Cij;
%     Alpha = params_local.NRTL.alpha;
% 
%     N = length(x);
%     Gamma = zeros(1,N);
%     Tau = Aij + Bij ./ T;
%     Gij = exp(-Alpha .* Tau);
% 
%     for i = 1:N
%         SUMC = 0; SUMD = 0; SUME = 0;
%         for j = 1:N
%             A = x(j) .* Gij(i,j);
%             SUMA = 0; SUMB = 0;
%             for k = 1:N
%                 SUMA = SUMA + x(k) .* Gij(k,j);
%                 SUMB = SUMB + x(k) .* Tau(k,j) .* Gij(k,j);
%             end
%             SUMA = max(SUMA,1e-16);
%             SUME = SUME + x(j) .* Gij(j,i);
%             SUMC = SUMC + A ./ SUMA .* (Tau(i,j) - SUMB ./ SUMA);
%             SUMD = SUMD + x(j) .* Tau(j,i) .* Gij(j,i);
%         end
%         SUME = max(SUME,1e-16);
%         Gamma(i) = exp(SUMD ./ SUME + SUMC);
%     end
% 
%     gamma = Gamma;
% end
