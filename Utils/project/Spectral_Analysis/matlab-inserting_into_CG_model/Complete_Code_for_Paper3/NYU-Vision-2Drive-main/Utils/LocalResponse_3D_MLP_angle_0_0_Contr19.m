%% LocalResponse_3D_MLP_angle_0_0_Contr19
% Routes to MLP regressors (angle_0_0_Contr19) produced from mlp_*_angle_0_0_Contr19.pt.
% LGNInd: [1,2,3,4] <-> [0,45,90,135] deg, 5 = background.
function LDEOutLIBy = LocalResponse_3D_MLP_angle_0_0_Contr19(celltype, L4EUse, L4IUse, LGNInd, L6LibInd)
    switch celltype
        case 'S'
            switch LGNInd
                case 1, LDEOutLIBy = S_LGNc1_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = S_LGNc2_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = S_LGNc3_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = S_LGNc4_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = S_LGNc5_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                otherwise, error('Invalid LGNInd for S');
            end
        case 'C'
            switch LGNInd
                case 1, LDEOutLIBy = C_LGNc1_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = C_LGNc2_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = C_LGNc3_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = C_LGNc4_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = C_LGNc5_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                otherwise, error('Invalid LGNInd for C');
            end
        case 'I'
            switch LGNInd
                case 1, LDEOutLIBy = I_LGNc1_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = I_LGNc2_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = I_LGNc3_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = I_LGNc4_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = I_LGNc5_MLP_angle_0_0_Contr19(L4EUse, L4IUse, L6LibInd);
                otherwise, error('Invalid LGNInd for I');
            end
        otherwise
            error('Invalid celltype');
    end
end
