% Analytic Jacobian of one-step map using 4D MLPs (angle_0_0, all contrasts, full dataset).
% Follows the same structure as compute_J_mlp_angle_0_0_phi_allContrasts but calls
% the *_MLP4D_angle_0_0_allContrasts_full_tanh_128_with_grad functions via
% local_mlp_grads_allContrasts_full_tanh_128.
%
% Inputs:
%   LDEUse        : struct with fields S, C, I (column vectors, length N)
%   ContrastUse   : vector of contrasts (same shape as LDEUse.S)
%   PixLGNCtgr    : N x 5 LGN mixing weights
%   L6Kernel      : convolution kernel used in L6 path
%   L6pars        : parameters for L6Convert (shared drive)
%   C_*_mean      : connectivity matrices (sparse)
%   L4*Ep         : normalizing scalars
%   N_HCOutY,NPixX,NPixY : spatial layout (rows = N_HCOutY*NPixY)
%   p             : (optional) convex mixing, x_{k+1}=p*F(x_k)+(1-p)*x_k
%   EKp, IKp      : (optional) L6Convert params for S/C saturation, and
%                   InhMulp params for I saturation
%   InhKillFlag   : (optional) enable the S/C/I saturation stage (default true)
%
% Output:
%   J : sparse Jacobian (3N x 3N)
function J = compute_J_mlp_angle_0_0_phi_allContrasts_full_tanh_128(LDEUse, ContrastUse, PixLGNCtgr, L6Kernel, L6pars, ...
    C_SS_mean, C_CS_mean, C_IS_mean, ...
    C_SC_mean, C_CC_mean, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, p, EKp, IKp, InhKillFlag)


    if nargin < 27 || isempty(InhKillFlag)
        InhKillFlag = true;
    end

    % ---------- Dimensions / helpers ----------
    wC = 0.3077;
    wS = 1 - wC;

    N = numel(LDEUse.S);
    Ny = N_HCOutY * NPixY;
    Nx = N / Ny;

    % Build circular convolution matrix for L6Kernel (persistent cache).
    K = build_conv_matrix_circular(L6Kernel, Ny, Nx);

    % ---------- InhKill preprocessing (elementwise) ----------
    E_raw = wS * LDEUse.S + wC * LDEUse.C;

    if InhKillFlag
        E_base = L6Convert(E_raw, EKp);
        E_adj = E_base ./ (E_raw);
        S_use = LDEUse.S .* E_adj;
        C_use = LDEUse.C .* E_adj;
        [I_use, ~] = InhMulp(LDEUse.I, IKp);
    else
        S_use = LDEUse.S;
        C_use = LDEUse.C;
        I_use = LDEUse.I;
    end

    E_use = wS * S_use + wC * C_use;

    L4EUse_S = (C_SS_mean * S_use + C_SC_mean * C_use) / L4SEp;
    L4EUse_C = (C_CS_mean * S_use + C_CC_mean * C_use) / L4CEp;
    L4EUse_I = (C_IS_mean * S_use + C_IC_mean * C_use) / L4IEp;

    L4IUse_S = (C_SI_mean * I_use) / L4SIp;
    L4IUse_C = (C_CI_mean * I_use) / L4CIp;
    L4IUse_I = (C_II_mean * I_use) / L4IIp;
            
    dL4E_S_S = (C_SS_mean) / L4SEp;
    dL4E_S_C = (C_SC_mean) / L4SEp;

    dL4E_C_S = (C_CS_mean) / L4CEp;
    dL4E_C_C = (C_CC_mean) / L4CEp;

    dL4E_I_S = (C_IS_mean) / L4IEp;
    dL4E_I_C = (C_IC_mean) / L4IEp;

    dL4I_S_I = (C_SI_mean) / L4SIp;
    dL4I_C_I = (C_CI_mean) / L4CIp;
    dL4I_I_I = (C_II_mean) / L4IIp;

    C_conv = K * E_use; % conv2 circular, vectorized
    L6EUse = L6Convert(reshape(C_conv, Ny, Nx), L6pars);
    L6EUse = L6EUse(:);

    L6grad = L6Convert_grad(reshape(C_conv, Ny, Nx), L6pars);
    L6grad = L6grad(:);

    lenMesh = max(ceil(max(L6EUse(:))/3), 2);
    L6ELibInd = L6EUse / 3;
    L6ELibInd(L6ELibInd < 1) = 1;
    L6ELibInd(L6ELibInd > lenMesh) = lenMesh;
    mask_L6 = (L6ELibInd > 1) & (L6ELibInd < lenMesh);

    DL6_base = spdiags((mask_L6/3) .* L6grad, 0, N, N);
    dL6_dS = DL6_base * K * (wS);
    dL6_dC = DL6_base * K * (wC);

    % ---------- Local MLP (per cell type, 4D all-contrast, full) ----------
    angle_tag = '0_0';
    [dPhiS_dL4E, dPhiS_dL4I, dPhiS_dL6] = local_mlp_grads_allContrasts_full_tanh_128('S', ...
        L4EUse_S, L4IUse_S, L6ELibInd, ContrastUse, PixLGNCtgr, angle_tag);
    [dPhiC_dL4E, dPhiC_dL4I, dPhiC_dL6] = local_mlp_grads_allContrasts_full_tanh_128('C', ...
        L4EUse_C, L4IUse_C, L6ELibInd, ContrastUse, PixLGNCtgr, angle_tag);
    [dPhiI_dL4E, dPhiI_dL4I, dPhiI_dL6] = local_mlp_grads_allContrasts_full_tanh_128('I', ...
        L4EUse_I, L4IUse_I, L6ELibInd, ContrastUse, PixLGNCtgr, angle_tag);

    DS_L4E = spdiags(dPhiS_dL4E, 0, N, N);
    DS_L4I = spdiags(dPhiS_dL4I, 0, N, N);
    DS_L6  = spdiags(dPhiS_dL6, 0, N, N);

    DC_L4E = spdiags(dPhiC_dL4E, 0, N, N);
    DC_L4I = spdiags(dPhiC_dL4I, 0, N, N);
    DC_L6  = spdiags(dPhiC_dL6, 0, N, N);

    DI_L4E = spdiags(dPhiI_dL4E, 0, N, N);
    DI_L4I = spdiags(dPhiI_dL4I, 0, N, N);
    DI_L6  = spdiags(dPhiI_dL6, 0, N, N);

    % ---------- Assemble Jacobian blocks ----------
    D_SS = DS_L4E * dL4E_S_S + DS_L4I * 0 + DS_L6 * dL6_dS;
    D_SC = DS_L4E * dL4E_S_C + DS_L4I * 0 + DS_L6 * dL6_dC;
    D_SI = DS_L4E * 0           + DS_L4I * dL4I_S_I + DS_L6 * 0;

    D_CS = DC_L4E * dL4E_C_S + DC_L4I * 0 + DC_L6 * dL6_dS;
    D_CC = DC_L4E * dL4E_C_C + DC_L4I * 0 + DC_L6 * dL6_dC;
    D_CI = DC_L4E * 0           + DC_L4I * dL4I_C_I + DC_L6 * 0;

    D_IS = DI_L4E * dL4E_I_S + DI_L4I * 0 + DI_L6 * dL6_dS;
    D_IC = DI_L4E * dL4E_I_C + DI_L4I * 0 + DI_L6 * dL6_dC;
    D_II = DI_L4E * 0           + DI_L4I * dL4I_I_I + DI_L6 * 0;

    J_core = [D_SS, D_SC, D_SI;
              D_CS, D_CC, D_CI;
              D_IS, D_IC, D_II];

  
        J = J_core;
    end


%% ---------- Helper: circular convolution matrix ----------
function K = build_conv_matrix_circular(kernel, Ny, Nx)
    persistent cache;
    key = sprintf('k%dx%d_%dx%d', Ny, Nx, size(kernel,1), size(kernel,2));
    if ~isempty(cache) && isfield(cache, key)
        K = cache.(key);
        return;
    end

    N = Ny * Nx;
    idx_grid = reshape(1:N, Ny, Nx);

    [kh, kw] = size(kernel);
    ctr_r = ceil(kh/2);
    ctr_c = ceil(kw/2);

    rows = [];
    cols = [];
    vals = [];

    for r = 1:kh
        for c = 1:kw
            weight = kernel(r,c);
            if weight == 0, continue; end
            dr = r - ctr_r;
            dc = c - ctr_c;
            src = circshift(idx_grid, [dr, dc]);
            rows = [rows; idx_grid(:)];
            cols = [cols; src(:)];
            vals = [vals; weight * ones(N,1)];
        end
    end

    K = sparse(rows, cols, vals, N, N);
    if isempty(cache)
        cache = struct();
    end
    cache.(key) = K;
end
