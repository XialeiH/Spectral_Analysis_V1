function [dPhi_dL4E, dPhi_dL4I, dPhi_dL6] = local_mlp_grads(celltype, L4EUse, L4IUse, L6Use, PixLGNCtgr, angle_tag)
    % Calls the appropriate *_MLP_angle_<angle>_with_grad for LGN 1-5 and accumulates with PixLGNCtgr.
    % angle_tag: string like '0_0', '7_5', '15_0', '22_5'
    N = length(L4EUse);
    dPhi_dL4E = zeros(N,1);
    dPhi_dL4I = zeros(N,1);
    dPhi_dL6 = zeros(N,1);
    for LGNInd = 1:5
        w = PixLGNCtgr(:,LGNInd);
        switch celltype
            case 'S'
                fn = str2func(sprintf('S_LGNc%d_MLP_angle_%s_with_grad', LGNInd, angle_tag));
            case 'C'
                fn = str2func(sprintf('C_LGNc%d_MLP_angle_%s_with_grad', LGNInd, angle_tag));
            case 'I'
                fn = str2func(sprintf('I_LGNc%d_MLP_angle_%s_with_grad', LGNInd, angle_tag));
            otherwise
                error('Invalid celltype');
        end
        [~, gE, gI, gL6] = fn(L4EUse, L4IUse, L6Use);
        dPhi_dL4E = dPhi_dL4E + w .* gE;
        dPhi_dL4I = dPhi_dL4I + w .* gI;
        dPhi_dL6 = dPhi_dL6 + w .* gL6;
    end
end
