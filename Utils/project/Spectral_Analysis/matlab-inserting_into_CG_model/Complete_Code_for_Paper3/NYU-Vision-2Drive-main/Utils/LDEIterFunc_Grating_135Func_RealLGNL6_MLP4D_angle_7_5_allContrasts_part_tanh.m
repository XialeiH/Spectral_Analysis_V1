%% 4D MLP angle_7_5 full version: uses LocalResponse_4D_MLP_angle_7_5_allContrasts_part_tanh
% Inputs must share shape; ContrastUse is the fourth input.
function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_MLP4D_angle_7_5_allContrasts_part_tanh(...
    celltype, L4EUse, L4IUse, ContrastUse, PixLGNCtgr, L6ELibInd)

    if ~isequal(size(L4EUse), size(L4IUse), size(ContrastUse))
        error('L4EUse, L4IUse, ContrastUse must have the same shape.');
    end

    LDEOutLIBy = cell(5,1);
    LibyAll = zeros(length(L4EUse),5); % 5 LGN inputs
    for LGNInd = 1:5
        LDEOutLIBy{LGNInd} = LocalResponse_4D_MLP_angle_7_5_allContrasts_part_tanh(...
            celltype, L4EUse, L4IUse, L6ELibInd, ContrastUse, LGNInd);
        LibyAll(:,LGNInd) = LDEOutLIBy{LGNInd};
    end
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);
end
