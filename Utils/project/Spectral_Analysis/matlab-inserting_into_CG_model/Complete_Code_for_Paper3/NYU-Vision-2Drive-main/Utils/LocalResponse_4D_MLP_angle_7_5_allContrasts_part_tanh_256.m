function LDEOutLIBy = LocalResponse_4D_MLP_angle_7_5_allContrasts_part_tanh_256(celltype, L4EUse, L4IUse, L6Use, ContrastUse, LGNInd) %#codegen
switch celltype
    case 'S'
        switch LGNInd
            case 1, LDEOutLIBy = S_LGNc1_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 2, LDEOutLIBy = S_LGNc2_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 3, LDEOutLIBy = S_LGNc3_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 4, LDEOutLIBy = S_LGNc4_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 5, LDEOutLIBy = S_LGNc5_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            otherwise, error('LGNInd out of range')
        end
    case 'C'
        switch LGNInd
            case 1, LDEOutLIBy = C_LGNc1_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 2, LDEOutLIBy = C_LGNc2_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 3, LDEOutLIBy = C_LGNc3_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 4, LDEOutLIBy = C_LGNc4_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 5, LDEOutLIBy = C_LGNc5_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            otherwise, error('LGNInd out of range')
        end
    case 'I'
        switch LGNInd
            case 1, LDEOutLIBy = I_LGNc1_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 2, LDEOutLIBy = I_LGNc2_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 3, LDEOutLIBy = I_LGNc3_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 4, LDEOutLIBy = I_LGNc4_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            case 5, LDEOutLIBy = I_LGNc5_MLP4D_angle_7_5_allContrasts_part_tanh_256(L4EUse, L4IUse, L6Use, ContrastUse);
            otherwise, error('LGNInd out of range')
        end
    otherwise, error('Unknown celltype')
end
end