function [D_SS,D_SC,D_SI, ...
          D_CS,D_CC,D_CI, ...
          D_IS,D_IC,D_II] = ...
    Compute_D_QR_from_2D_sigmoid_chain_rule(f, ParaODE)


    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
    N         = PixNumOut;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

    S = f(idxS);
    C = f(idxC);
    I = f(idxI);

    % connectivity
    C_SS_mean = ParaODE.C_SS_mean;
    C_CS_mean = ParaODE.C_CS_mean;
    C_IS_mean = ParaODE.C_IS_mean;

    C_SC_mean = ParaODE.C_SC_mean;
    C_CC_mean = ParaODE.C_CC_mean;
    C_IC_mean = ParaODE.C_IC_mean;

    C_SI_mean = ParaODE.C_SI_mean;
    C_CI_mean = ParaODE.C_CI_mean;
    C_II_mean = ParaODE.C_II_mean;

    PixInptCtgrUse = ParaODE.PixInptCtgrUse;  % N×4×4


    L4EUse = C_SS_mean*S + ...
             C_CS_mean*S + ...
             C_IS_mean*S + ...
             C_SC_mean*C + ...
             C_CC_mean*C + ...
             C_IC_mean*C;

    L4IUse = C_SI_mean*I + ...
             C_CI_mean*I + ...
             C_II_mean*I;

    L4EUse = L4EUse(:);
    L4IUse = L4IUse(:);


    h = 1e-4;

    % constant for each block
    c_SS = C_SS_mean + C_CS_mean + C_IS_mean;
    c_CS = C_SS_mean + C_CS_mean + C_IS_mean;
    c_IS = C_SS_mean + C_CS_mean + C_IS_mean;

    c_SC = C_SC_mean + C_CC_mean + C_IC_mean;
    c_CC = C_SC_mean + C_CC_mean + C_IC_mean;
    c_IC = C_SC_mean + C_CC_mean + C_IC_mean;

    c_SI = C_SI_mean + C_CI_mean + C_II_mean;
    c_CI = C_SI_mean + C_CI_mean + C_II_mean;
    c_II = C_SI_mean + C_CI_mean + C_II_mean;


%%
    D_SS = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_S_plus  = local_Psi('S', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_S_minus = local_Psi('S', L4E_minus, L4IUse, PixInptCtgrUse);

        D_SS(:,j) = (Psi_S_plus - Psi_S_minus) / (2*h);
        
    end
    D_SS = D_SS * c_SS;

%%
    D_SC = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_S_plus  = local_Psi('S', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_S_minus = local_Psi('S', L4E_minus, L4IUse, PixInptCtgrUse);

        D_SC(:,j) = (Psi_S_plus - Psi_S_minus) / (2*h);
        
    end
    D_SC = D_SC * c_SC;

%%
    D_SI = zeros(N, N);
    for j = 1:N
        L4I_plus  = L4IUse;
        L4I_minus = L4IUse;
        L4I_plus(j)  = L4IUse(j) + h;
        L4I_minus(j) = L4IUse(j) - h;

        Psi_S_plus  = local_Psi('S', L4EUse,  L4I_plus,  PixInptCtgrUse);
        Psi_S_minus = local_Psi('S', L4EUse, L4I_minus, PixInptCtgrUse);

        D_SI(:,j) = (Psi_S_plus - Psi_S_minus) / (2*h);
        
    end
    D_SI = D_SI * c_SI;

%%
    D_CS = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_C_plus  = local_Psi('C', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_C_minus = local_Psi('C', L4E_minus, L4IUse, PixInptCtgrUse);

        D_CS(:,j) = (Psi_C_plus - Psi_C_minus) / (2*h);
        
    end
    D_CS = D_CS * c_CS;

%%
    D_CC = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_C_plus  = local_Psi('C', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_C_minus = local_Psi('C', L4E_minus, L4IUse, PixInptCtgrUse);

        D_CC(:,j) = (Psi_C_plus - Psi_C_minus) / (2*h);
        
    end
    D_CC = D_CC * c_CC;

%%
    D_CI = zeros(N, N);
    for j = 1:N
        L4I_plus  = L4IUse;
        L4I_minus = L4IUse;
        L4I_plus(j)  = L4IUse(j) + h;
        L4I_minus(j) = L4IUse(j) - h;

        Psi_C_plus  = local_Psi('C', L4EUse,  L4I_plus,  PixInptCtgrUse);
        Psi_C_minus = local_Psi('C', L4EUse, L4I_minus, PixInptCtgrUse);

        D_CI(:,j) = (Psi_C_plus - Psi_C_minus) / (2*h);
        
    end
    D_CI = D_CI * c_CI;

%%
    D_IS = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_I_plus  = local_Psi('I', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_I_minus = local_Psi('I', L4E_minus, L4IUse, PixInptCtgrUse);

        D_IS(:,j) = (Psi_I_plus - Psi_I_minus) / (2*h);
        
    end
    D_IS = D_IS * c_IS;

%%
    D_IC = zeros(N, N);
    for j = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4E_plus(j)  = L4EUse(j) + h;
        L4E_minus(j) = L4EUse(j) - h;

        Psi_I_plus  = local_Psi('I', L4E_plus,  L4IUse,  PixInptCtgrUse);
        Psi_I_minus = local_Psi('I', L4E_minus, L4IUse, PixInptCtgrUse);

        D_IC(:,j) = (Psi_I_plus - Psi_I_minus) / (2*h);
        
    end
    D_IC = D_IC * c_IC;

%%
    D_II = zeros(N, N);
    for j = 1:N
        L4I_plus  = L4IUse;
        L4I_minus = L4IUse;
        L4I_plus(j)  = L4IUse(j) + h;
        L4I_minus(j) = L4IUse(j) - h;

        Psi_I_plus  = local_Psi('I', L4EUse,  L4I_plus,  PixInptCtgrUse);
        Psi_I_minus = local_Psi('I', L4EUse, L4I_minus, PixInptCtgrUse);

        D_II(:,j) = (Psi_I_plus - Psi_I_minus) / (2*h);
        
    end
    D_II = D_II * c_II;

end



function Psi_Q = local_Psi(celltype, L4EUse, L4IUse, PixInptCtgrUse)

    Npix = numel(L4EUse);
    LibyAll = zeros(Npix,4,4);  % [pixel, LGN, L6]

    switch celltype
        case 'S'
            % LGN c1
            LibyAll(:,1,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,2) = S_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,3) = S_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,4) = S_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
            % LGN c2
            LibyAll(:,2,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,2) = S_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,3) = S_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,4) = S_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
            % LGN c3
            LibyAll(:,3,1) = S_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,2) = S_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,3) = S_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,4) = S_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
            % LGN c4
            LibyAll(:,4,1) = S_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,2) = S_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,3) = S_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,4) = S_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);

        case 'C'
            LibyAll(:,1,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,2) = C_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,3) = C_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,4) = C_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,2,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,2) = C_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,3) = C_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,4) = C_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,3,1) = C_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,2) = C_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,3) = C_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,4) = C_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,4,1) = C_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,2) = C_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,3) = C_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,4) = C_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);

        case 'I'
            LibyAll(:,1,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,2) = I_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,3) = I_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,4) = I_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,2,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,2) = I_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,3) = I_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,4) = I_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,3,1) = I_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,2) = I_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,3) = I_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,4) = I_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);

            LibyAll(:,4,1) = I_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,2) = I_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,3) = I_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,4,4) = I_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
    end

    % Pixelwise weighted sum over the 16 LGN–L6 combinations.
    % PixInptCtgrUse is N×4×4, LibyAll is N×4×4.
    Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2,3]);
    Psi_Q = Psi_Q(:);  % ensure column vector N×1
end
