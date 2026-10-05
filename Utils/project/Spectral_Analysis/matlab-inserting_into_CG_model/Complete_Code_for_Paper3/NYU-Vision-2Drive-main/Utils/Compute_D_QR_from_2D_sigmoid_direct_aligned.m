function D = Compute_D_QR_from_2D_sigmoid_direct_aligned(f, ParaODE)
% Finite-difference Jacobian using the direct sigmoid path, aligned to the
% iteration-library settings (same h and pixelwise construction only).
%
% Differences vs Compute_D_QR_from_2D_sigmoid_direct:
%   - step size h = 1e-3 (matches iteration-library FD)
%   - otherwise identical one-shot forward using *_Sigmoid.

    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
    N         = PixNumOut;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

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

    h = 1e-3;  % align with iteration-library FD step

    D = zeros(3*N, 3*N);

    for j = 1:3*N
        f_plus = f;
        f_minus = f;

        f_plus(j)  = f_plus(j)  + h;
        f_minus(j) = f_minus(j) - h;

        S_plus = f_plus(idxS); C_plus = f_plus(idxC); I_plus = f_plus(idxI);
        S_minus = f_minus(idxS); C_minus = f_minus(idxC); I_minus = f_minus(idxI);

        L4E_plus = C_SS_mean*S_plus + C_CS_mean*S_plus + C_IS_mean*S_plus + ...
                   C_SC_mean*C_plus + C_CC_mean*C_plus + C_IC_mean*C_plus;
        L4I_plus = C_SI_mean*I_plus + C_CI_mean*I_plus + C_II_mean*I_plus;

        L4E_minus = C_SS_mean*S_minus + C_CS_mean*S_minus + C_IS_mean*S_minus + ...
                    C_SC_mean*C_minus + C_CC_mean*C_minus + C_IC_mean*C_minus;
        L4I_minus = C_SI_mean*I_minus + C_CI_mean*I_minus + C_II_mean*I_minus;

        Psi_S_p = local_Psi('S', L4E_plus,  L4I_plus,  PixInptCtgrUse);
        Psi_S_m = local_Psi('S', L4E_minus, L4I_minus, PixInptCtgrUse);

        Psi_C_p = local_Psi('C', L4E_plus,  L4I_plus,  PixInptCtgrUse);
        Psi_C_m = local_Psi('C', L4E_minus, L4I_minus, PixInptCtgrUse);

        Psi_I_p = local_Psi('I', L4E_plus,  L4I_plus,  PixInptCtgrUse);
        Psi_I_m = local_Psi('I', L4E_minus, L4I_minus, PixInptCtgrUse);

        D(:, j) = [Psi_S_p - Psi_S_m;
                   Psi_C_p - Psi_C_m;
                   Psi_I_p - Psi_I_m] / (2*h);
    end
end


function Psi_Q = local_Psi(celltype, L4EUse, L4IUse, PixInptCtgrUse)
    Npix = numel(L4EUse);
    LibyAll = zeros(Npix,4,4);  % [pixel, LGN, L6]

    switch celltype
        case 'S'
            LibyAll(:,1,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,2) = S_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,3) = S_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,1,4) = S_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,2) = S_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,3) = S_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,2,4) = S_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,1) = S_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,2) = S_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,3) = S_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
            LibyAll(:,3,4) = S_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
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

    Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2,3]);
    Psi_Q = Psi_Q(:);
end
