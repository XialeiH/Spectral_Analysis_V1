function LDEOutLIBy = LocalResponse_5D_MLP5D_allOrientaion_allContrasts_part_tanh_relu(celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, LGNInd) %#codegen
switch celltype
    case 'S'
        switch LGNInd
            case 1, LDEOutLIBy = S_LGNc1_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 2, LDEOutLIBy = S_LGNc2_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 3, LDEOutLIBy = S_LGNc3_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 4, LDEOutLIBy = S_LGNc4_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 5, LDEOutLIBy = S_LGNc5_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            otherwise, error('LGNInd out of range')
        end
    case 'C'
        switch LGNInd
            case 1, LDEOutLIBy = C_LGNc1_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 2, LDEOutLIBy = C_LGNc2_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 3, LDEOutLIBy = C_LGNc3_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 4, LDEOutLIBy = C_LGNc4_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 5, LDEOutLIBy = C_LGNc5_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            otherwise, error('LGNInd out of range')
        end
    case 'I'
        switch LGNInd
            case 1, LDEOutLIBy = I_LGNc1_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 2, LDEOutLIBy = I_LGNc2_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 3, LDEOutLIBy = I_LGNc3_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 4, LDEOutLIBy = I_LGNc4_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            case 5, LDEOutLIBy = I_LGNc5_MLP5D_allOrientaion_allContrasts_part_tanh_relu(L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse);
            otherwise, error('LGNInd out of range')
        end
    otherwise, error('Unknown celltype')
end
end
