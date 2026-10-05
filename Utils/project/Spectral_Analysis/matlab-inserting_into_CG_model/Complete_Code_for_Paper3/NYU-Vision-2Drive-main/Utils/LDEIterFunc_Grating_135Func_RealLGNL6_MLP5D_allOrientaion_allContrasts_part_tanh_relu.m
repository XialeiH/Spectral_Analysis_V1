%% 5D MLP full version: uses LocalResponse_5D_MLP5D_allOrientaion_allContrasts_part_tanh_relu
function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_MLP5D_allOrientaion_allContrasts_part_tanh_relu(...
    celltype, L4EUse, L4IUse, ContrastUse, OrientationUse, PixLGNCtgr, L6ELibInd)

    if ~isequal(size(L4EUse), size(L4IUse), size(ContrastUse), size(OrientationUse))
        error('L4EUse, L4IUse, ContrastUse, OrientationUse must have the same shape.');
    end

    LDEOutLIBy = cell(5,1);
    LibyAll = zeros(length(L4EUse),5); %% 5 LGN inputs
    for LGNInd = 1:5
        LDEOutLIBy{LGNInd} = LocalResponse_5D_MLP5D_allOrientaion_allContrasts_part_tanh_relu(...
            celltype, L4EUse, L4IUse, L6ELibInd, ContrastUse, OrientationUse, LGNInd);
        LibyAll(:,LGNInd) = LDEOutLIBy{LGNInd};
    end
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);
end
