function [DpPhi, paramNames, detail] = GA_h96_pref6D_parameter_dphi(LDEUse, ContrastUse, OrientationUse, PixLGNCtgr, L6Kernel, L6pars, ...
    C_SS_mean, C_CS_mean, C_IS_mean, ...
    C_SC_mean, C_CC_mean, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, EKp, IKp, InhKillFlag)
% Direct derivatives DpPhi for local log-strength circuit parameters.
%
% The parameter coordinates are log gains. For example, d/dp for
% S_from_C is the response change caused by C_SC_mean -> exp(p)*C_SC_mean
% at p = 0. The NN itself is fixed; parameters only change the L4E/L4I/L6
% inputs seen by the local response function.
if nargin < 26 || isempty(InhKillFlag)
    InhKillFlag = true;
end

spec = GA_h96_parameter_spec();
paramNames = spec.Names;

wC = 0.3077;
wS = 1 - wC;
N = numel(LDEUse.S);
ContrastUse = expand_to_state(ContrastUse, N, 'ContrastUse');
OrientationUse = expand_to_state(OrientationUse, N, 'OrientationUse');

if InhKillFlag
    E_raw = wS * LDEUse.S + wC * LDEUse.C;
    E_base = L6Convert(E_raw, EKp);
    E_adj = E_base ./ E_raw;
    E_adj(~isfinite(E_adj)) = 0;
    S_use = LDEUse.S .* E_adj;
    C_use = LDEUse.C .* E_adj;
    I_use = InhMulp(LDEUse.I, IKp);
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

FieldRow = N_HCOutY * NPixY;
FieldCol = floor(N / FieldRow);
L4Efield = reshape(E_use, FieldRow, FieldCol);
L4Efield_padded = padarray(L4Efield, [1, 1], 'circular');
Cconv = conv2(L4Efield_padded, L6Kernel, 'same');
Cconv = Cconv(2:end-1, 2:end-1);
L6EUse = L6Convert(Cconv, L6pars);

lenMesh = 40;
L6ELibIndRaw = L6EUse(:) / 3;
L6ELibInd = L6ELibIndRaw;
L6ELibInd(L6ELibInd < 1) = 1;
L6ELibInd(L6ELibInd > lenMesh) = lenMesh;
maskL6 = (L6ELibIndRaw > 1) & (L6ELibIndRaw < lenMesh);
L6grad = L6Convert_grad(Cconv, L6pars);
dL6Lib_dLogGain = (maskL6 / 3) .* L6grad(:) .* Cconv(:);

[dPhiS_dL4E, dPhiS_dL4I, dPhiS_dL6] = local_h96baseline_pref6D_grads('S', L4EUse_S, L4IUse_S, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);
[dPhiC_dL4E, dPhiC_dL4I, dPhiC_dL6] = local_h96baseline_pref6D_grads('C', L4EUse_C, L4IUse_C, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);
[dPhiI_dL4E, dPhiI_dL4I, dPhiI_dL6] = local_h96baseline_pref6D_grads('I', L4EUse_I, L4IUse_I, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);

DpPhi = zeros(3*N, spec.NumParameters);
Sidx = 1:N;
Cidx = N + (1:N);
Iidx = 2*N + (1:N);

% E-to-E paths.
DpPhi(Sidx, 1) = dPhiS_dL4E .* ((C_SS_mean * S_use) / L4SEp);
DpPhi(Sidx, 2) = dPhiS_dL4E .* ((C_SC_mean * C_use) / L4SEp);
DpPhi(Cidx, 3) = dPhiC_dL4E .* ((C_CS_mean * S_use) / L4CEp);
DpPhi(Cidx, 4) = dPhiC_dL4E .* ((C_CC_mean * C_use) / L4CEp);

% E-to-I paths.
DpPhi(Iidx, 5) = dPhiI_dL4E .* ((C_IS_mean * S_use) / L4IEp);
DpPhi(Iidx, 6) = dPhiI_dL4E .* ((C_IC_mean * C_use) / L4IEp);

% I-to-E and I-to-I paths.
DpPhi(Sidx, 7) = dPhiS_dL4I .* ((C_SI_mean * I_use) / L4SIp);
DpPhi(Cidx, 8) = dPhiC_dL4I .* ((C_CI_mean * I_use) / L4CIp);
DpPhi(Iidx, 9) = dPhiI_dL4I .* ((C_II_mean * I_use) / L4IIp);

% Common L6 feedback input gain.
DpPhi(Sidx, 10) = dPhiS_dL6 .* dL6Lib_dLogGain;
DpPhi(Cidx, 10) = dPhiC_dL6 .* dL6Lib_dLogGain;
DpPhi(Iidx, 10) = dPhiI_dL6 .* dL6Lib_dLogGain;

detail = struct();
detail.ParameterSpec = spec;
detail.L6Use = L6ELibInd;
detail.L6RawIndex = L6ELibIndRaw;
detail.L6InputDerivative = dL6Lib_dLogGain;
detail.ContrastUse = ContrastUse;
detail.OrientationUse = OrientationUse;
detail.Columns = paramNames;
end

function x = expand_to_state(x, n, name)
if isscalar(x)
    x = x * ones(n, 1);
else
    x = x(:);
end
if numel(x) ~= n
    error('%s must be scalar or have %d elements.', name, n);
end
end
