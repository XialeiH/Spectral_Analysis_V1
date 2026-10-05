function [Psi_Q, dPsi_dL4E, dPsi_dL4I] = local_Psi_with_deriv(celltype, L4EUse, L4IUse, PixInptCtgrUse)
%LOCAL_PSI_WITH_DERIV Weighted LGN->L6 response and derivatives.
%   Mirrors local_Psi.m but calls the *_with_deriv sigmoid variants so the
%   weighted output Psi_Q and its sensitivities to L4E/L4I can be obtained
%   from the same PixInptCtgrUse tensor.

    Npix = numel(L4EUse);
    LibyAll = zeros(Npix, 4, 4);
    dLiby_dL4E = zeros(Npix, 4, 4);
    dLiby_dL4I = zeros(Npix, 4, 4);

    funcNames = getSigmoidFunctionNames(celltype);

    for c = 1:4
        for l = 1:4
            funcHandle = str2func(funcNames{c,l});
            [LibyAll(:,c,l), dLiby_dL4E(:,c,l), dLiby_dL4I(:,c,l)] = ...
                funcHandle(L4EUse, L4IUse);
        end
    end

    Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2, 3]);
    Psi_Q = Psi_Q(:);

    dPsi_dL4E = sum(PixInptCtgrUse .* dLiby_dL4E, [2, 3]);
    dPsi_dL4E = dPsi_dL4E(:);

    dPsi_dL4I = sum(PixInptCtgrUse .* dLiby_dL4I, [2, 3]);
    dPsi_dL4I = dPsi_dL4I(:);
end


function funcNames = getSigmoidFunctionNames(celltype)
%GETSIGMOIDFUNCTIONNAMES Return the 4x4 table of sigmoid function names.

    switch celltype
        case 'S'
            funcNames = {
                'S_LGNc1_L6c1_Sigmoid_with_deriv', 'S_LGNc1_L6c2_Sigmoid_with_deriv', 'S_LGNc1_L6c3_Sigmoid_with_deriv', 'S_LGNc1_L6c4_Sigmoid_with_deriv';
                'S_LGNc1_L6c1_Sigmoid_with_deriv', 'S_LGNc2_L6c2_Sigmoid_with_deriv', 'S_LGNc2_L6c3_Sigmoid_with_deriv', 'S_LGNc2_L6c4_Sigmoid_with_deriv';
                'S_LGNc3_L6c1_Sigmoid_with_deriv', 'S_LGNc3_L6c2_Sigmoid_with_deriv', 'S_LGNc3_L6c3_Sigmoid_with_deriv', 'S_LGNc3_L6c4_Sigmoid_with_deriv';
                'S_LGNc4_L6c1_Sigmoid_with_deriv', 'S_LGNc4_L6c2_Sigmoid_with_deriv', 'S_LGNc4_L6c3_Sigmoid_with_deriv', 'S_LGNc4_L6c4_Sigmoid_with_deriv'};

        case 'C'
            funcNames = {
                'C_LGNc1_L6c1_Sigmoid_with_deriv', 'C_LGNc1_L6c2_Sigmoid_with_deriv', 'C_LGNc1_L6c3_Sigmoid_with_deriv', 'C_LGNc1_L6c4_Sigmoid_with_deriv';
                'C_LGNc1_L6c1_Sigmoid_with_deriv', 'C_LGNc2_L6c2_Sigmoid_with_deriv', 'C_LGNc2_L6c3_Sigmoid_with_deriv', 'C_LGNc2_L6c4_Sigmoid_with_deriv';
                'C_LGNc3_L6c1_Sigmoid_with_deriv', 'C_LGNc3_L6c2_Sigmoid_with_deriv', 'C_LGNc3_L6c3_Sigmoid_with_deriv', 'C_LGNc3_L6c4_Sigmoid_with_deriv';
                'C_LGNc4_L6c1_Sigmoid_with_deriv', 'C_LGNc4_L6c2_Sigmoid_with_deriv', 'C_LGNc4_L6c3_Sigmoid_with_deriv', 'C_LGNc4_L6c4_Sigmoid_with_deriv'};

        case 'I'
            funcNames = {
                'I_LGNc1_L6c1_Sigmoid_with_deriv', 'I_LGNc1_L6c2_Sigmoid_with_deriv', 'I_LGNc1_L6c3_Sigmoid_with_deriv', 'I_LGNc1_L6c4_Sigmoid_with_deriv';
                'I_LGNc1_L6c1_Sigmoid_with_deriv', 'I_LGNc2_L6c2_Sigmoid_with_deriv', 'I_LGNc2_L6c3_Sigmoid_with_deriv', 'I_LGNc2_L6c4_Sigmoid_with_deriv';
                'I_LGNc3_L6c1_Sigmoid_with_deriv', 'I_LGNc3_L6c2_Sigmoid_with_deriv', 'I_LGNc3_L6c3_Sigmoid_with_deriv', 'I_LGNc3_L6c4_Sigmoid_with_deriv';
                'I_LGNc4_L6c1_Sigmoid_with_deriv', 'I_LGNc4_L6c2_Sigmoid_with_deriv', 'I_LGNc4_L6c3_Sigmoid_with_deriv', 'I_LGNc4_L6c4_Sigmoid_with_deriv'};

        otherwise
            error('local_Psi_with_deriv:UnknownCellType', ...
                  'Cell type must be ''S'', ''C'', or ''I''.');
    end
end
