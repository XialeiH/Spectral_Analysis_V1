function [Psi_Q, dPsi_dL4EProd, dPsi_dL4IProd] = local_Psi_with_deriv_scaled(celltype, L4EUse, L4IUse, PixInptCtgrUse, L4QEp, L4QIp)
%LOCAL_PSI_WITH_DERIV_SCALED Weighted LGN->L6 response and derivatives
% w.r.t. products (L4EUse*L4QEp, L4IUse*L4QIp), where L4QEp/L4QIp are scalars.

    Npix = numel(L4EUse);
    LibyAll = zeros(Npix, 4, 4);
    dLiby_dL4EProd = zeros(Npix, 4, 4);
    dLiby_dL4IProd = zeros(Npix, 4, 4);

    funcNames = getSigmoidFunctionNames(celltype);

    for c = 1:4
        for l = 1:4
            funcHandle = str2func(funcNames{c,l});
            [LibyAll(:,c,l), dLiby_dL4EProd(:,c,l), dLiby_dL4IProd(:,c,l)] = ...
                funcHandle(L4EUse, L4IUse, L4QEp, L4QIp);
        end
    end

    Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2, 3]);
    Psi_Q = Psi_Q(:);

    dPsi_dL4EProd = sum(PixInptCtgrUse .* dLiby_dL4EProd, [2, 3]);
    dPsi_dL4EProd = dPsi_dL4EProd(:);

    dPsi_dL4IProd = sum(PixInptCtgrUse .* dLiby_dL4IProd, [2, 3]);
    dPsi_dL4IProd = dPsi_dL4IProd(:);
end


function funcNames = getSigmoidFunctionNames(celltype)
    switch celltype
        case 'S'
            funcNames = {
                'S_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'S_LGNc1_L6c2_Sigmoid_with_deriv_scaled', 'S_LGNc1_L6c3_Sigmoid_with_deriv_scaled', 'S_LGNc1_L6c4_Sigmoid_with_deriv_scaled';
                'S_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'S_LGNc2_L6c2_Sigmoid_with_deriv_scaled', 'S_LGNc2_L6c3_Sigmoid_with_deriv_scaled', 'S_LGNc2_L6c4_Sigmoid_with_deriv_scaled';
                'S_LGNc3_L6c1_Sigmoid_with_deriv_scaled', 'S_LGNc3_L6c2_Sigmoid_with_deriv_scaled', 'S_LGNc3_L6c3_Sigmoid_with_deriv_scaled', 'S_LGNc3_L6c4_Sigmoid_with_deriv_scaled';
                'S_LGNc4_L6c1_Sigmoid_with_deriv_scaled', 'S_LGNc4_L6c2_Sigmoid_with_deriv_scaled', 'S_LGNc4_L6c3_Sigmoid_with_deriv_scaled', 'S_LGNc4_L6c4_Sigmoid_with_deriv_scaled'};

        case 'C'
            funcNames = {
                'C_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'C_LGNc1_L6c2_Sigmoid_with_deriv_scaled', 'C_LGNc1_L6c3_Sigmoid_with_deriv_scaled', 'C_LGNc1_L6c4_Sigmoid_with_deriv_scaled';
                'C_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'C_LGNc2_L6c2_Sigmoid_with_deriv_scaled', 'C_LGNc2_L6c3_Sigmoid_with_deriv_scaled', 'C_LGNc2_L6c4_Sigmoid_with_deriv_scaled';
                'C_LGNc3_L6c1_Sigmoid_with_deriv_scaled', 'C_LGNc3_L6c2_Sigmoid_with_deriv_scaled', 'C_LGNc3_L6c3_Sigmoid_with_deriv_scaled', 'C_LGNc3_L6c4_Sigmoid_with_deriv_scaled';
                'C_LGNc4_L6c1_Sigmoid_with_deriv_scaled', 'C_LGNc4_L6c2_Sigmoid_with_deriv_scaled', 'C_LGNc4_L6c3_Sigmoid_with_deriv_scaled', 'C_LGNc4_L6c4_Sigmoid_with_deriv_scaled'};

        case 'I'
            funcNames = {
                'I_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'I_LGNc1_L6c2_Sigmoid_with_deriv_scaled', 'I_LGNc1_L6c3_Sigmoid_with_deriv_scaled', 'I_LGNc1_L6c4_Sigmoid_with_deriv_scaled';
                'I_LGNc1_L6c1_Sigmoid_with_deriv_scaled', 'I_LGNc2_L6c2_Sigmoid_with_deriv_scaled', 'I_LGNc2_L6c3_Sigmoid_with_deriv_scaled', 'I_LGNc2_L6c4_Sigmoid_with_deriv_scaled';
                'I_LGNc3_L6c1_Sigmoid_with_deriv_scaled', 'I_LGNc3_L6c2_Sigmoid_with_deriv_scaled', 'I_LGNc3_L6c3_Sigmoid_with_deriv_scaled', 'I_LGNc3_L6c4_Sigmoid_with_deriv_scaled';
                'I_LGNc4_L6c1_Sigmoid_with_deriv_scaled', 'I_LGNc4_L6c2_Sigmoid_with_deriv_scaled', 'I_LGNc4_L6c3_Sigmoid_with_deriv_scaled', 'I_LGNc4_L6c4_Sigmoid_with_deriv_scaled'};

        otherwise
            error('local_Psi_with_deriv_scaled:UnknownCellType', ...
                  'Cell type must be ''S'', ''C'', or ''I''.');
    end
end
