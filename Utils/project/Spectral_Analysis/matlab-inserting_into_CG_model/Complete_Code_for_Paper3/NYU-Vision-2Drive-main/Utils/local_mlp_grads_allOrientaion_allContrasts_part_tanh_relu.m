function [dPhi_dL4E, dPhi_dL4I, dPhi_dL6] = local_mlp_grads_allOrientaion_allContrasts_part_tanh_relu(celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, PixLGNCtgr)
% Accumulates analytical gradients for 5D MLP models.
LGNnum = size(PixLGNCtgr,2);
dPhi_dL4E = zeros(size(L4EUse));
dPhi_dL4I = zeros(size(L4IUse));
dPhi_dL6  = zeros(size(L6Use));
for LGNInd = 1:LGNnum
    switch celltype
        case 'S'
            fn = str2func(sprintf('S_LGNc%d_MLP5D_allOrientaion_allContrasts_part_tanh_relu_with_grad', LGNInd));
        case 'C'
            fn = str2func(sprintf('C_LGNc%d_MLP5D_allOrientaion_allContrasts_part_tanh_relu_with_grad', LGNInd));
        case 'I'
            fn = str2func(sprintf('I_LGNc%d_MLP5D_allOrientaion_allContrasts_part_tanh_relu_with_grad', LGNInd));
        otherwise
            error('Unknown celltype');
    end
    [~, gE, gI, gL6] = fn(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
    dPhi_dL4E = dPhi_dL4E + PixLGNCtgr(:,LGNInd).*gE;
    dPhi_dL4I = dPhi_dL4I + PixLGNCtgr(:,LGNInd).*gI;
    dPhi_dL6  = dPhi_dL6  + PixLGNCtgr(:,LGNInd).*gL6;
end
end
