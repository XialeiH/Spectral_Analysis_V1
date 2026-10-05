function [J, diagnostics] = compute_J_library135_phi(LDEUse, PixLGNCtgr, ...
    L6Kernel, L6pars, C_SS_mean, C_CS_mean, C_IS_mean, ...
    C_SC_mean, C_CC_mean, C_IC_mean, C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, CtgrOrderReadout, ...
    N_HCOutY, NPixX, NPixY, EKp, IKp, InhKillFlag)
% Jacobian of the genuine 135-function Paper 3 table-interpolation map.

% The linear network chain is differentiated analytically. Local response
% derivatives are evaluated from the same trilinearly interpolated tables
% used by LDEIteration_135FuncMain_CombDom_RealLGNL6.

% This is the derivative of F(x), not of the relaxed iteration
% p*F(x)+(1-p)*x. That matches the J_F convention used by the DNN row.

if nargin < 29 || isempty(InhKillFlag)
    InhKillFlag = true;
end

wC = 0.3077;
wS = 1 - wC;
N = numel(LDEUse.S);
Ny = N_HCOutY * NPixY;
Nx = N / Ny;
assert(Nx == round(Nx), 'The map dimensions are inconsistent.');
Nx = round(Nx);

K = local_build_conv_matrix_circular(L6Kernel, Ny, Nx);
Eraw = wS * LDEUse.S + wC * LDEUse.C;

if InhKillFlag
    Ebase = L6Convert(Eraw, EKp);
    Eadj = Ebase ./ Eraw;
    dEbase = L6Convert_grad(Eraw, EKp);
    dEadj = (dEbase .* Eraw - Ebase) ./ (Eraw .^ 2);
    dEadj(~isfinite(dEadj)) = 0;

    Suse = LDEUse.S .* Eadj;
    Cuse = LDEUse.C .* Eadj;
    Iuse = InhMulp(LDEUse.I, IKp);
    dIuse = local_inhmulp_grad(LDEUse.I, IKp);
else
    Eadj = ones(N, 1);
    dEadj = zeros(N, 1);
    Suse = LDEUse.S;
    Cuse = LDEUse.C;
    Iuse = LDEUse.I;
    dIuse = ones(N, 1);
end

D_Suse_S = spdiags(Eadj + LDEUse.S .* dEadj * wS, 0, N, N);
D_Suse_C = spdiags(LDEUse.S .* dEadj * wC, 0, N, N);
D_Cuse_S = spdiags(LDEUse.C .* dEadj * wS, 0, N, N);
D_Cuse_C = spdiags(Eadj + LDEUse.C .* dEadj * wC, 0, N, N);
D_Iuse_I = spdiags(dIuse, 0, N, N);

L4E_S = (C_SS_mean * Suse + C_SC_mean * Cuse) / L4SEp;
L4E_C = (C_CS_mean * Suse + C_CC_mean * Cuse) / L4CEp;
L4E_I = (C_IS_mean * Suse + C_IC_mean * Cuse) / L4IEp;
L4I_S = (C_SI_mean * Iuse) / L4SIp;
L4I_C = (C_CI_mean * Iuse) / L4CIp;
L4I_I = (C_II_mean * Iuse) / L4IIp;

dL4E_S_S = (C_SS_mean * D_Suse_S + C_SC_mean * D_Cuse_S) / L4SEp;
dL4E_S_C = (C_SS_mean * D_Suse_C + C_SC_mean * D_Cuse_C) / L4SEp;
dL4E_C_S = (C_CS_mean * D_Suse_S + C_CC_mean * D_Cuse_S) / L4CEp;
dL4E_C_C = (C_CS_mean * D_Suse_C + C_CC_mean * D_Cuse_C) / L4CEp;
dL4E_I_S = (C_IS_mean * D_Suse_S + C_IC_mean * D_Cuse_S) / L4IEp;
dL4E_I_C = (C_IS_mean * D_Suse_C + C_IC_mean * D_Cuse_C) / L4IEp;
dL4I_S_I = (C_SI_mean * D_Iuse_I) / L4SIp;
dL4I_C_I = (C_CI_mean * D_Iuse_I) / L4CIp;
dL4I_I_I = (C_II_mean * D_Iuse_I) / L4IIp;

Euse = wS * Suse + wC * Cuse;
dEuse_S = wS * D_Suse_S + wC * D_Cuse_S;
dEuse_C = wS * D_Suse_C + wC * D_Cuse_C;
Cconv = K * Euse;
L6Hz = L6Convert(reshape(Cconv, Ny, Nx), L6pars);
L6Hz = L6Hz(:);
L6indexRaw = L6Hz / 3;
L6meshCount = size(LDEFrfuncAll{1}.S, 2);
L6index = min(max(L6indexRaw, 1), L6meshCount);
L6curveGrad = L6Convert_grad(reshape(Cconv, Ny, Nx), L6pars);
L6curveGrad = L6curveGrad(:);
insideL6 = L6indexRaw > 1 & L6indexRaw < L6meshCount;
DL6 = spdiags((insideL6 / 3) .* L6curveGrad, 0, N, N);
dL6_S = DL6 * K * dEuse_S;
dL6_C = DL6 * K * dEuse_C;

[gSE, gSI, gS6, domainS] = local_table_grads('S', L4E_S, L4I_S, ...
    L6index, PixLGNCtgr, L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, CtgrOrderReadout);
[gCE, gCI, gC6, domainC] = local_table_grads('C', L4E_C, L4I_C, ...
    L6index, PixLGNCtgr, L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, CtgrOrderReadout);
[gIE, gII, gI6, domainI] = local_table_grads('I', L4E_I, L4I_I, ...
    L6index, PixLGNCtgr, L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, CtgrOrderReadout);

GSE = spdiags(gSE, 0, N, N); GSI = spdiags(gSI, 0, N, N); GS6 = spdiags(gS6, 0, N, N);
GCE = spdiags(gCE, 0, N, N); GCI = spdiags(gCI, 0, N, N); GC6 = spdiags(gC6, 0, N, N);
GIE = spdiags(gIE, 0, N, N); GII = spdiags(gII, 0, N, N); GI6 = spdiags(gI6, 0, N, N);

J = [GSE*dL4E_S_S + GS6*dL6_S, GSE*dL4E_S_C + GS6*dL6_C, GSI*dL4I_S_I; ...
     GCE*dL4E_C_S + GC6*dL6_S, GCE*dL4E_C_C + GC6*dL6_C, GCI*dL4I_C_I; ...
     GIE*dL4E_I_S + GI6*dL6_S, GIE*dL4E_I_C + GI6*dL6_C, GII*dL4I_I_I];

diagnostics = struct('DomainS', domainS, 'DomainC', domainC, ...
    'DomainI', domainI, 'L6IndexMin', min(L6index), ...
    'L6IndexMax', max(L6index), 'L6ClampCount', nnz(~insideL6));
end

function [dE, dI, dL6, domainIndex] = local_table_grads(celltype, ...
    xq, yq, zq, weights, xAll, yAll, tableAll, categoryOrder)
% Evaluate central derivatives of MATLAB's piecewise-linear interp3 map.
% The steps are tiny fractions of each table spacing, so the result is the
% exact local slope except at a table knot, where it is the symmetric slope.

for domainIndex = 1:numel(xAll)
    xGrid = unique(xAll{domainIndex}(:));
    yGrid = unique(yAll{domainIndex}(:));
    zGrid = (1:size(tableAll{domainIndex}.S, 2))';
    tableMatrix = local_pack_tables(tableAll{domainIndex}, celltype, size(weights, 2));

    hx = max(min(diff(xGrid)) * 1e-4, 1e-7);
    hy = max(min(diff(yGrid)) * 1e-4, 1e-7);
    hz = max(min(diff(zGrid)) * 1e-4, 1e-7);
    dE = (local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq + hx, yq, zq, weights, categoryOrder) - ...
        local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq - hx, yq, zq, weights, categoryOrder)) / (2 * hx);
    dI = (local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq, yq + hy, zq, weights, categoryOrder) - ...
        local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq, yq - hy, zq, weights, categoryOrder)) / (2 * hy);
    dL6 = (local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq, yq, zq + hz, weights, categoryOrder) - ...
        local_table_response(xGrid, yGrid, zGrid, tableMatrix, ...
        xq, yq, zq - hz, weights, categoryOrder)) / (2 * hz);

    if all(isfinite([dE; dI; dL6]))
        return;
    end
end
error('No library domain contains every equilibrium query and perturbation.');
end

function matrix = local_pack_tables(data, celltype, categoryCount)
source = data.(celltype);
l6Count = size(source, 2);
sample = source{1, 1};
matrix = zeros(categoryCount, size(sample, 1), size(sample, 2), l6Count, 'like', sample);
for category = 1:categoryCount
    for l6 = 1:l6Count
        matrix(category, :, :, l6) = source{category, l6};
    end
end
end

function response = local_table_response(xGrid, yGrid, zGrid, matrix, ...
    xq, yq, zq, weights, categoryOrder)
[XX, YY, ZZ] = meshgrid(xGrid, yGrid, zGrid);
categoryCount = size(weights, 2);
values = zeros(numel(xq), categoryCount);
for category = 1:categoryCount
    values(:, category) = interp3(XX, YY, ZZ, ...
        squeeze(matrix(category, :, :, :)), xq, yq, zq, 'linear');
end
if categoryCount >= 4
    foreground = values(:, 1:4);
    values(:, 1:4) = foreground(:, categoryOrder(:)');
end
response = sum(weights .* values, 2);
end

function derivative = local_inhmulp_grad(x, parameters)
h = 1e-5 * max(1, abs(x));
derivative = (InhMulp(x + h, parameters) - InhMulp(x - h, parameters)) ./ (2 * h);
derivative(~isfinite(derivative)) = 0;
end

function K = local_build_conv_matrix_circular(kernel, rowsCount, columnsCount)
N = rowsCount * columnsCount;
indices = reshape(1:N, rowsCount, columnsCount);
[kernelRows, kernelColumns] = size(kernel);
centerRow = ceil(kernelRows / 2);
centerColumn = ceil(kernelColumns / 2);
rows = [];
columns = [];
values = [];
for row = 1:kernelRows
    for column = 1:kernelColumns
        weight = kernel(row, column);
        if weight == 0
            continue;
        end
        source = circshift(indices, [row - centerRow, column - centerColumn]);
        rows = [rows; indices(:)]; %#ok<AGROW>
        columns = [columns; source(:)]; %#ok<AGROW>
        values = [values; weight * ones(N, 1)]; %#ok<AGROW>
    end
end
K = sparse(rows, columns, values, N, N);
end
