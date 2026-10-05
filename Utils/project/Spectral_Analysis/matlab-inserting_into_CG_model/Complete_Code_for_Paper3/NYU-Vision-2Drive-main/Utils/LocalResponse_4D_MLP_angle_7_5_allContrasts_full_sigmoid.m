%% LocalResponse_4D_MLP_angle_7_5_allContrasts_full_sigmoid
% Routes to 4D MLP regressors (angle_7_5, all contrasts, full dataset).
% Inputs must share the same shape; ContrastUse is the fourth input (scaled by 1/100 inside the MLP call).
% LGNInd: [1,2,3,4] <-> [0,45,90,135] deg, 5 = background.
function LDEOutLIBy = LocalResponse_4D_MLP_angle_7_5_allContrasts_full_sigmoid(celltype, L4EUse, L4IUse, L6Use, ContrastUse, LGNInd)
    if ~isequal(size(L4EUse), size(L4IUse), size(L6Use), size(ContrastUse))
        error('Inputs must have the same shape.');
    end

    switch celltype
        case 'S'
            switch LGNInd
                case 1, LDEOutLIBy = S_LGNc1_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 2, LDEOutLIBy = S_LGNc2_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 3, LDEOutLIBy = S_LGNc3_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 4, LDEOutLIBy = S_LGNc4_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 5, LDEOutLIBy = S_LGNc5_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                otherwise, error('Invalid LGNInd for S celltype');
            end

        case 'C'
            switch LGNInd
                case 1, LDEOutLIBy = C_LGNc1_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 2, LDEOutLIBy = C_LGNc2_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 3, LDEOutLIBy = C_LGNc3_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 4, LDEOutLIBy = C_LGNc4_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 5, LDEOutLIBy = C_LGNc5_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                otherwise, error('Invalid LGNInd for C celltype');
            end

        case 'I'
            switch LGNInd
                case 1, LDEOutLIBy = I_LGNc1_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 2, LDEOutLIBy = I_LGNc2_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 3, LDEOutLIBy = I_LGNc3_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 4, LDEOutLIBy = I_LGNc4_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                case 5, LDEOutLIBy = I_LGNc5_MLP4D_angle_7_5_allContrasts_full_sigmoid(L4EUse, L4IUse, L6Use, ContrastUse);
                otherwise, error('Invalid LGNInd for I celltype');
            end
        otherwise
            error('Invalid celltype');
    end
end
