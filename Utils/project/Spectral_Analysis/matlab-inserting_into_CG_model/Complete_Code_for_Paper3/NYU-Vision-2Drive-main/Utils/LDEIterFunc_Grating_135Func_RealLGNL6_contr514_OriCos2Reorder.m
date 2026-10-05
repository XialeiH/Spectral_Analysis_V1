%% Iterations of LDE: Use precomputed function to determine output of a cell type
% Variant: for contrast 5/14, reorder canonical library LGN1..4 outputs by
% OriCos2 computed from the preferred angle and the input angle.
%
% Canonical foreground LGN categories use preferred angles:
%   LGN1 -> 0 deg
%   LGN2 -> 135 deg
%   LGN3 -> 90 deg
%   LGN4 -> 45 deg
%
% For contrast 5/14 only, we compute OriCos2 for these four categories and
% sort the foreground channels in descending OriCos2 order before combining
% them with PixLGNCtgr. Background channel 5 is left unchanged.

function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_contr514_OriCos2Reorder(...
    L4EmeshX, L4ImeshY, L6MeshZ, ...
    LDEFrfunc_Subf, ...
    L4EUse, L4IUse, ...
    PixLGNCtgr, L6ELibInd, ...
    AlphaUse, ContrastUse)

    LDEOutLIBy = cell(size(LDEFrfunc_Subf, 1), 1);
    LibyAll = zeros(length(L4EUse), size(LDEFrfunc_Subf, 1));
    [XX, YY, ZZ] = meshgrid(unique(L4EmeshX), unique(L4ImeshY), L6MeshZ);
    for LGNInd = 1:size(LDEFrfunc_Subf, 1)
        LDEOutLIBy{LGNInd} = interp3(XX, YY, ZZ, ...
            squeeze(LDEFrfunc_Subf(LGNInd, :, :, :)), ...
            L4EUse, L4IUse, L6ELibInd, 'linear');
        LibyAll(:, LGNInd) = LDEOutLIBy{LGNInd};
    end

    contrastScalar = local_scalar_input(ContrastUse, 'ContrastUse');
    angleScalar = local_scalar_input(AlphaUse, 'AlphaUse');
    if ismember(round(contrastScalar), [5, 14])
        thetaMap = [0, 135, 90, 45];
        gammaUse = abs(local_wrap_to90(angleScalar - thetaMap));
        oriCos2 = cosd(2 .* gammaUse);
        % Reorder canonical LGN1..4 by descending OriCos2; break ties by
        % canonical category id so the behavior is deterministic.
        sortTable = [(-oriCos2(:)), (1:4)'];
        sortTable = sortrows(sortTable, [1, 2]);
        reorderIdx = sortTable(:, 2)';
        LibyFG = LibyAll(:, 1:4);
        LibyAll(:, 1:4) = LibyFG(:, reorderIdx);
    end

    LDEOutS = sum(PixLGNCtgr .* LibyAll, 2);
end

function scalarValue = local_scalar_input(x, name)
    x = x(:);
    if isempty(x)
        error('%s must not be empty.', name);
    end
    scalarValue = double(x(1));
    if any(abs(double(x) - scalarValue) > 1e-9)
        error('%s must be scalar or spatially constant for the OriCos2 reorder variant.', name);
    end
end

function thetaWrapped = local_wrap_to90(thetaDeg)
    thetaWrapped = mod(thetaDeg + 90, 180) - 90;
end
