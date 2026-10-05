%% LocalResponse_3D_moe_v5
% Routes to v5 MoE sigmoids (deep gate + normalized residual tanh features).
% LGNInd: [1,2,3,4] <-> [0,45,90,135] deg, 5 = background.
function LDEOutLIBy = LocalResponse_3D_moe_v5(celltype, L4EUse, L4IUse, LGNInd, L6LibInd)
    switch celltype
        case 'S'
            switch LGNInd
                case 1, LDEOutLIBy = S_LGNc1_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = S_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = S_LGNc3_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = S_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = S_LGNc5_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
            end

        case 'C'
            switch LGNInd
                case 1, LDEOutLIBy = C_LGNc1_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = C_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = C_LGNc3_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = C_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = C_LGNc5_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
            end
        case 'I'
            switch LGNInd
                case 1, LDEOutLIBy = I_LGNc1_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 2, LDEOutLIBy = I_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 3, LDEOutLIBy = I_LGNc3_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 4, LDEOutLIBy = I_LGNc2_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
                case 5, LDEOutLIBy = I_LGNc5_3DSigmoid_moe_v5(L4EUse, L4IUse, L6LibInd);
            end
    end
end
