%% MoE version: uses LocalResponse_3D_moe (soft mixture sigmoids) instead of the hard-segment sigmoids
% Zhuo-Cheng Xiao style wrapper, adapted to call *_3DSigmoid_moe functions.

function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_Sigmoid_moe(...
    celltype, L4EUse, L4IUse, PixLGNCtgr, L6ELibInd)

    LDEOutLIBy = cell(5,1);
    LibyAll = zeros(length(L4EUse),5); % 5 LGN inputs
    
    % LGN1-5
    LDEOutLIBy{1} = LocalResponse_3D_moe(celltype, L4EUse, L4IUse, 1, L6ELibInd);
    LibyAll(:,1) = LDEOutLIBy{1};
   
    LDEOutLIBy{2} = LocalResponse_3D_moe(celltype, L4EUse, L4IUse, 2, L6ELibInd);
    LibyAll(:,2) = LDEOutLIBy{2};

    LDEOutLIBy{3} = LocalResponse_3D_moe(celltype, L4EUse, L4IUse, 3, L6ELibInd);
    LibyAll(:,3) = LDEOutLIBy{3};

    LDEOutLIBy{4} = LocalResponse_3D_moe(celltype, L4EUse, L4IUse, 4, L6ELibInd);
    LibyAll(:,4) = LDEOutLIBy{4};

    LDEOutLIBy{5} = LocalResponse_3D_moe(celltype, L4EUse, L4IUse, 5, L6ELibInd);
    LibyAll(:,5) = LDEOutLIBy{5};
   
    % Weighted combine over LGN inputs
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);
end
