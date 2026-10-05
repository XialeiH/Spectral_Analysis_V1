function J = compute_J_h96baseline_pref6D_phi_perturbed_3range(LDEUse, ContrastUse, OrientationUse, PixLGNCtgr, L6Kernel, L6pars, ...
    C_SS_mean, C_CS_mean, C_IS_mean, ...
    C_SC_mean, C_CC_mean, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, p, EKp, IKp, InhKillFlag, L6EquBlendWeight, Perturb)
if nargin < 28 || isempty(InhKillFlag); InhKillFlag = true; end
if nargin < 29 || isempty(L6EquBlendWeight); L6EquBlendWeight = 0; end
if nargin < 30 || isempty(Perturb); Perturb = struct(); end
L6EquBlendWeight = double(L6EquBlendWeight);
S_EL6_scale = getfield_default(Perturb, 'S_EL6', 1);
S_IL6_scale = getfield_default(Perturb, 'S_IL6', 1);
L4InputParam = getfield_default(Perturb, 'L4InputParam', '');
L4InputScale = getfield_default(Perturb, 'L4InputScale', 1);
L4InputScaleDeriv = getfield_default(Perturb, 'L4InputScaleDeriv', true);
wC = 0.3077; wS = 1 - wC;
N = numel(LDEUse.S);
Ny = round(double(N_HCOutY) * double(NPixY));
Nx_guess = round(double(NPixX));
if Ny <= 0
    error('compute_J_h96baseline_pref6D_phi:InvalidNy', 'Computed Ny must be positive. Got Ny=%g.', Ny);
end
if Ny * Nx_guess == N
    Nx = Nx_guess;
else
    Nx_from_N = round(double(N) / double(Ny));
    if Ny * Nx_from_N ~= N
        error('compute_J_h96baseline_pref6D_phi:InconsistentGrid', ...
            'Cannot infer integer Ny/Nx: N=%d, Ny=%d, NPixX=%g, NPixY=%g, N_HCOutY=%g.', ...
            N, Ny, double(NPixX), double(NPixY), double(N_HCOutY));
    end
    Nx = Nx_from_N;
end
K = build_conv_matrix_circular(L6Kernel, Ny, Nx);
E_raw = wS * LDEUse.S + wC * LDEUse.C;
if InhKillFlag
    E_base = L6Convert(E_raw, EKp);
    E_adj = E_base ./ (E_raw);
    dEbase_dEraw = L6Convert_grad(E_raw, EKp);
    dEadj_dEraw = (dEbase_dEraw .* E_raw - E_base) ./ (E_raw.^2);
    dEadj_dEraw(~isfinite(dEadj_dEraw)) = 0;
    S_use = LDEUse.S .* E_adj; C_use = LDEUse.C .* E_adj; [I_use, ~] = InhMulp(LDEUse.I, IKp);
    dIuse_dI = InhMulp_grad(LDEUse.I, IKp);
else
    S_use = LDEUse.S; C_use = LDEUse.C; I_use = LDEUse.I;
    dEadj_dEraw = zeros(N,1);
    E_adj = ones(N,1);
    dIuse_dI = ones(N,1);
end
E_use = wS * S_use + wC * C_use;
L4EUse_S = (C_SS_mean * S_use + C_SC_mean * C_use) / L4SEp;
L4EUse_C = (C_CS_mean * S_use + C_CC_mean * C_use) / L4CEp;
L4EUse_I = (C_IS_mean * S_use + C_IC_mean * C_use) / L4IEp;
L4IUse_S = (C_SI_mean * I_use) / L4SIp;
L4IUse_C = (C_CI_mean * I_use) / L4CIp;
L4IUse_I = (C_II_mean * I_use) / L4IIp;
D_Suse_S = spdiags(E_adj + LDEUse.S .* dEadj_dEraw * wS, 0, N, N);
D_Suse_C = spdiags(LDEUse.S .* dEadj_dEraw * wC, 0, N, N);
D_Cuse_S = spdiags(LDEUse.C .* dEadj_dEraw * wS, 0, N, N);
D_Cuse_C = spdiags(E_adj + LDEUse.C .* dEadj_dEraw * wC, 0, N, N);
D_Iuse_I = spdiags(dIuse_dI, 0, N, N);
dL4E_S_S = (C_SS_mean * D_Suse_S + C_SC_mean * D_Cuse_S) / L4SEp;
dL4E_S_C = (C_SS_mean * D_Suse_C + C_SC_mean * D_Cuse_C) / L4SEp;
dL4E_C_S = (C_CS_mean * D_Suse_S + C_CC_mean * D_Cuse_S) / L4CEp;
dL4E_C_C = (C_CS_mean * D_Suse_C + C_CC_mean * D_Cuse_C) / L4CEp;
dL4E_I_S = (C_IS_mean * D_Suse_S + C_IC_mean * D_Cuse_S) / L4IEp;
dL4E_I_C = (C_IS_mean * D_Suse_C + C_IC_mean * D_Cuse_C) / L4IEp;
dL4I_S_I = (C_SI_mean * D_Iuse_I) / L4SIp;
dL4I_C_I = (C_CI_mean * D_Iuse_I) / L4CIp;
dL4I_I_I = (C_II_mean * D_Iuse_I) / L4IIp;
if ~isempty(L4InputParam) && L4InputScale ~= 1
    switch L4InputParam
        case 'S_EE'
            L4EUse_S = L4EUse_S * L4InputScale;
            L4EUse_C = L4EUse_C * L4InputScale;
            if L4InputScaleDeriv
                dL4E_S_S = dL4E_S_S * L4InputScale;
                dL4E_S_C = dL4E_S_C * L4InputScale;
                dL4E_C_S = dL4E_C_S * L4InputScale;
                dL4E_C_C = dL4E_C_C * L4InputScale;
            end
        case 'S_IE'
            L4EUse_I = L4EUse_I * L4InputScale;
            if L4InputScaleDeriv
                dL4E_I_S = dL4E_I_S * L4InputScale;
                dL4E_I_C = dL4E_I_C * L4InputScale;
            end
        case 'S_EI'
            L4IUse_S = L4IUse_S * L4InputScale;
            L4IUse_C = L4IUse_C * L4InputScale;
            if L4InputScaleDeriv
                dL4I_S_I = dL4I_S_I * L4InputScale;
                dL4I_C_I = dL4I_C_I * L4InputScale;
            end
        case 'S_II'
            L4IUse_I = L4IUse_I * L4InputScale;
            if L4InputScaleDeriv
                dL4I_I_I = dL4I_I_I * L4InputScale;
            end
        otherwise
            error('compute_J_h96baseline_pref6D_phi_perturbed:UnknownL4InputParam', ...
                'Unknown L4InputParam %s.', L4InputParam);
    end
end
C_conv = K * E_use; L6EUse = L6Convert(reshape(C_conv, Ny, Nx), L6pars); L6EUse = L6EUse(:);
lenMesh = 40;
L6ELibIndRaw_E = S_EL6_scale * (L6EUse / 3);
L6ELibIndRaw_I = S_IL6_scale * (L6EUse / 3);
L6ELibInd_E = L6ELibIndRaw_E; L6ELibInd_E(L6ELibInd_E < 1) = 1; L6ELibInd_E(L6ELibInd_E > lenMesh) = lenMesh;
L6ELibInd_I = L6ELibIndRaw_I; L6ELibInd_I(L6ELibInd_I < 1) = 1; L6ELibInd_I(L6ELibInd_I > lenMesh) = lenMesh;
mask_L6_E = (L6ELibIndRaw_E > 1) & (L6ELibIndRaw_E < lenMesh);
mask_L6_I = (L6ELibIndRaw_I > 1) & (L6ELibIndRaw_I < lenMesh);
L6grad = L6Convert_grad(reshape(C_conv, Ny, Nx), L6pars); L6grad = L6grad(:);
DL6_base_E = spdiags((mask_L6_E * S_EL6_scale / 3) .* L6grad, 0, N, N);
DL6_base_I = spdiags((mask_L6_I * S_IL6_scale / 3) .* L6grad, 0, N, N);
dEuse_dS = wS * D_Suse_S + wC * D_Cuse_S;
dEuse_dC = wS * D_Suse_C + wC * D_Cuse_C;
% New L6 Jacobian knob:
%   L6 = w*L6equ + (1-w)*L6original.
% At the fixed point, L6equ equals the current L6 value but is fixed with
% respect to state perturbations, so only the original L6 path contributes
% derivative weight (1-w).
dL6E_dS = (1 - L6EquBlendWeight) * (DL6_base_E * K * dEuse_dS);
dL6E_dC = (1 - L6EquBlendWeight) * (DL6_base_E * K * dEuse_dC);
dL6I_dS = (1 - L6EquBlendWeight) * (DL6_base_I * K * dEuse_dS);
dL6I_dC = (1 - L6EquBlendWeight) * (DL6_base_I * K * dEuse_dC);
[dPhiS_dL4E, dPhiS_dL4I, dPhiS_dL6] = local_h96baseline_pref6D_grads('S', L4EUse_S, L4IUse_S, L6ELibInd_E, ContrastUse, OrientationUse, PixLGNCtgr);
[dPhiC_dL4E, dPhiC_dL4I, dPhiC_dL6] = local_h96baseline_pref6D_grads('C', L4EUse_C, L4IUse_C, L6ELibInd_E, ContrastUse, OrientationUse, PixLGNCtgr);
if getfield_default(Perturb, 'MomentCorrection', false)
    [dPhiI_dL4E, dPhiI_dL4I, dPhiI_dL6] = local_h96baseline_pref6D_grads_momentcorrected( ...
        'I', L4EUse_I, L4IUse_I, L6ELibInd_I, ...
        ContrastUse, OrientationUse, PixLGNCtgr, Perturb);
else
    [dPhiI_dL4E, dPhiI_dL4I, dPhiI_dL6] = local_h96baseline_pref6D_grads( ...
        'I', L4EUse_I, L4IUse_I, L6ELibInd_I, ...
        ContrastUse, OrientationUse, PixLGNCtgr);
end
DS_L4E = spdiags(dPhiS_dL4E, 0, N, N); DS_L4I = spdiags(dPhiS_dL4I, 0, N, N); DS_L6 = spdiags(dPhiS_dL6, 0, N, N);
DC_L4E = spdiags(dPhiC_dL4E, 0, N, N); DC_L4I = spdiags(dPhiC_dL4I, 0, N, N); DC_L6 = spdiags(dPhiC_dL6, 0, N, N);
DI_L4E = spdiags(dPhiI_dL4E, 0, N, N); DI_L4I = spdiags(dPhiI_dL4I, 0, N, N); DI_L6 = spdiags(dPhiI_dL6, 0, N, N);
D_SS = DS_L4E * dL4E_S_S + DS_L4I * 0 + DS_L6 * dL6E_dS;
D_SC = DS_L4E * dL4E_S_C + DS_L4I * 0 + DS_L6 * dL6E_dC;
D_SI = DS_L4E * 0           + DS_L4I * dL4I_S_I + DS_L6 * 0;
D_CS = DC_L4E * dL4E_C_S + DC_L4I * 0 + DC_L6 * dL6E_dS;
D_CC = DC_L4E * dL4E_C_C + DC_L4I * 0 + DC_L6 * dL6E_dC;
D_CI = DC_L4E * 0           + DC_L4I * dL4I_C_I + DC_L6 * 0;
D_IS = DI_L4E * dL4E_I_S + DI_L4I * 0 + DI_L6 * dL6I_dS;
D_IC = DI_L4E * dL4E_I_C + DI_L4I * 0 + DI_L6 * dL6I_dC;
D_II = DI_L4E * 0           + DI_L4I * dL4I_I_I + DI_L6 * 0;
J_core = [D_SS, D_SC, D_SI; D_CS, D_CC, D_CI; D_IS, D_IC, D_II];

    J = J_core;
end

function K = build_conv_matrix_circular(kernel, Ny, Nx)
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
end

function val = getfield_default(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end

function dy = InhMulp_grad(x, IKp)
h = 1e-5 * max(1, abs(x));
yPlus = InhMulp(x + h, IKp);
yMinus = InhMulp(x - h, IKp);
dy = (yPlus - yMinus) ./ (2*h);
dy(~isfinite(dy)) = 0;
end
