function [dPhi_dL4E, dPhi_dL4I, dPhi_dL6] = local_h96baseline_pref6D_grads(celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, PixLGNCtgr)
LGNnum = size(PixLGNCtgr, 2);
dPhi_dL4E = zeros(size(L4EUse));
dPhi_dL4I = zeros(size(L4IUse));
dPhi_dL6 = zeros(size(L6Use));

for LGNInd = 1:LGNnum
    LGNctgrUse = LGNInd * ones(size(OrientationUse));
    [~, gE, gI, gL6] = LocalResponse_6D_MLP_prefAngle_with_grad(celltype, OrientationUse, LGNctgrUse, L4EUse, L4IUse, L6Use, ContrastUse);
    dPhi_dL4E = dPhi_dL4E + PixLGNCtgr(:, LGNInd) .* gE;
    dPhi_dL4I = dPhi_dL4I + PixLGNCtgr(:, LGNInd) .* gI;
    dPhi_dL6 = dPhi_dL6 + PixLGNCtgr(:, LGNInd) .* gL6;
end
end
