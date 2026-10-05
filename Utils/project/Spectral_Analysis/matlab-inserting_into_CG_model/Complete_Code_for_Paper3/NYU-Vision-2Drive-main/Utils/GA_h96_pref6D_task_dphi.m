function [dPhi, featureNames, detail] = GA_h96_pref6D_task_dphi(LDEUse, ContrastUse, OrientationUse, PixLGNCtgr, L6Kernel, L6pars, ...
    C_SS_mean, C_CS_mean, C_IS_mean, ...
    C_SC_mean, C_CC_mean, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, EKp, IKp, InhKillFlag)
% Direct task derivatives dPhi/du for the h96 pref6D NN at a fixed state.
% Columns are derivatives with respect to [log_contrast, orientation_deg].
if nargin < 26 || isempty(InhKillFlag)
    InhKillFlag = true;
end

wC = 0.3077;
wS = 1 - wC;
N = numel(LDEUse.S);
if NPixX <= 0 || NPixY <= 0
    error('NPixX and NPixY must be positive.');
end
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
L6EUse = L6Convert(Cconv(2:end-1, 2:end-1), L6pars);

lenMesh = 40;
L6ELibInd = L6EUse(:) / 3;
L6ELibInd(L6ELibInd < 1) = 1;
L6ELibInd(L6ELibInd > lenMesh) = lenMesh;

[dS_dLogC, dS_dOri] = local_task_grads('S', L4EUse_S, L4IUse_S, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);
[dC_dLogC, dC_dOri] = local_task_grads('C', L4EUse_C, L4IUse_C, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);
[dI_dLogC, dI_dOri] = local_task_grads('I', L4EUse_I, L4IUse_I, L6ELibInd, ContrastUse, OrientationUse, PixLGNCtgr);

dPhi = [dS_dLogC, dS_dOri; ...
        dC_dLogC, dC_dOri; ...
        dI_dLogC, dI_dOri];
featureNames = {'log_contrast', 'orientation_deg'};

detail = struct();
detail.L6Use = L6ELibInd;
detail.ContrastUse = ContrastUse;
detail.OrientationUse = OrientationUse;
detail.Columns = featureNames;
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

function [dLogContrast, dOrientation] = local_task_grads(celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, PixLGNCtgr)
n = numel(L4EUse);
dContrast = zeros(n, 1);
dOrientation = zeros(n, 1);
thetaMap = [0, 135, 90, 45];

for LGNInd = 1:size(PixLGNCtgr, 2)
    fgMask = (LGNInd >= 1) && (LGNInd <= 4);
    contrastEff = ContrastUse;
    oriCos2 = ones(n, 1);
    dOriCos2_dTheta = zeros(n, 1);

    if fgMask
        thetaWrapped = wrapTo90_pref(OrientationUse - thetaMap(LGNInd));
        oriCos2 = cosd(2 .* thetaWrapped);
        dOriCos2_dTheta = -(pi / 90) .* sind(2 .* thetaWrapped);
    else
        contrastEff(:) = 0;
    end

    switch celltype
        case 'S'
            [~, ~, ~, ~, gContrast, gOriCos2] = predict_pref6D_S_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
        case 'C'
            [~, ~, ~, ~, gContrast, gOriCos2] = predict_pref6D_C_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
        case 'I'
            [~, ~, ~, ~, gContrast, gOriCos2] = predict_pref6D_I_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
        otherwise
            error('Unknown celltype %s.', celltype);
    end

    w = PixLGNCtgr(:, LGNInd);
    dContrast = dContrast + w .* gContrast .* double(fgMask);
    dOrientation = dOrientation + w .* gOriCos2 .* dOriCos2_dTheta;
end

dLogContrast = ContrastUse .* dContrast;
end

function thetaWrapped = wrapTo90_pref(thetaDeg)
thetaWrapped = mod(thetaDeg + 90, 180) - 90;
end
