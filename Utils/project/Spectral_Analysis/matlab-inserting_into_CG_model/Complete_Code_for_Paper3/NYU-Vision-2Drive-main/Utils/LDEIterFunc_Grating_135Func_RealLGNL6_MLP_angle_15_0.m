%% MLP angle_15_0 version: uses LocalResponse_3D_MLP_angle_15_0 (tanh-MLP regressors)
function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_MLP_angle_15_0(...
    celltype, L4EUse, L4IUse, PixLGNCtgr, L6ELibInd)

    LDEOutLIBy = cell(5,1);
    LibyAll = zeros(length(L4EUse),5); % 5 LGN inputs
    for LGNInd = 1:5
        LDEOutLIBy{LGNInd} = LocalResponse_3D_MLP_angle_15_0(celltype, L4EUse, L4IUse, LGNInd, L6ELibInd);
        LibyAll(:,LGNInd) = LDEOutLIBy{LGNInd};
    end
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);
end
