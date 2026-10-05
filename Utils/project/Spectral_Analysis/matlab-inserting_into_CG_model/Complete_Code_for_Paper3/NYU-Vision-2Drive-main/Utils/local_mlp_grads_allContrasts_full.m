function [dPhi_dL4E, dPhi_dL4I, dPhi_dL6] = local_mlp_grads_allContrasts_full(celltype, L4EUse, L4IUse, L6Use, ContrastUse, PixLGNCtgr, angle_tag)
    % Calls the appropriate *_MLP4D_angle_<angle>_allContrasts_full_with_grad for LGN 1-5 and accumulates with PixLGNCtgr.
    N = length(L4EUse);
    dPhi_dL4E = zeros(N,1);
    dPhi_dL4I = zeros(N,1);
    dPhi_dL6 = zeros(N,1);
    for LGNInd = 1:5
        w = PixLGNCtgr(:,LGNInd);
        switch celltype
            case 'S'
                fn = str2func(sprintf('S_LGNc%d_MLP4D_angle_%s_allContrasts_full_with_grad', LGNInd, angle_tag));
            case 'C'
                fn = str2func(sprintf('C_LGNc%d_MLP4D_angle_%s_allContrasts_full_with_grad', LGNInd, angle_tag));
            case 'I'
                fn = str2func(sprintf('I_LGNc%d_MLP4D_angle_%s_allContrasts_full_with_grad', LGNInd, angle_tag));
            otherwise
                error('Invalid celltype');
        end
        [~, gE, gI, gL6] = fn(L4EUse, L4IUse, L6Use, ContrastUse);
        dPhi_dL4E = dPhi_dL4E + w .* gE;
        dPhi_dL4I = dPhi_dL4I + w .* gI;
        dPhi_dL6 = dPhi_dL6 + w .* gL6;
    end
end
