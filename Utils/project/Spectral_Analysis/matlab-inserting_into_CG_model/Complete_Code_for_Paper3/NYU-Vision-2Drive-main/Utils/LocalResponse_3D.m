%% this function is for knowing local response of each of the 3200 pixels, calling corresponding Sigmoid functions
% for now only consider the input fr in the small domain, the NN is trained
% by data within the small domain. 
% % LGNInd is based on the OS in each HC 
% % LGNInd: [1,2,3,4] <-> [0,45,90,135] degree, LGN 5: background LGN input


function LDEOutLIBy = LocalResponse_3D(celltype, L4EUse, L4IUse, LGNInd, L6LibInd)
    switch celltype
        case 'S'
            switch LGNInd
                case 1, LDEOutLIBy = S_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 2, LDEOutLIBy = S_LGNc2_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 3, LDEOutLIBy = S_LGNc3_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 4, LDEOutLIBy = S_LGNc2_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 5, LDEOutLIBy = S_LGNc5_3DSigmoid(L4EUse, L4IUse,L6LibInd);
            end

        case 'C'
            switch LGNInd
                case 1, LDEOutLIBy = C_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 2, LDEOutLIBy = C_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 3, LDEOutLIBy = C_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 4, LDEOutLIBy = C_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 5, LDEOutLIBy = C_LGNc5_3DSigmoid(L4EUse, L4IUse,L6LibInd);
            end
        case 'I'
            switch LGNInd
                case 1, LDEOutLIBy = I_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 2, LDEOutLIBy = I_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 3, LDEOutLIBy = I_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 4, LDEOutLIBy = I_LGNc1_3DSigmoid(L4EUse, L4IUse,L6LibInd);
                case 5, LDEOutLIBy = I_LGNc5_3DSigmoid(L4EUse, L4IUse,L6LibInd);
            end
    end
end