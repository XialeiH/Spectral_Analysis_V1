function D = Compute_D_QR_from_2D_sigmoid_direct(f, ParaODE)
% Compute_D_QR_from_2D_sigmoid
%
% For each postsynaptic population Q ∈ {S,C,I} and each presynaptic
% population R ∈ {S,C,I}, compute
%
%   D_QR(i,j) = ∂Ψ_Q(i) / ∂(R-input from pixel j)
%
% using finite differences in the *presynaptic* activity at pixel j and
% recomputing L4E/L4I and the local response Ψ_Q for all pixels.

    %% Basic sizes and indexing
    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
    N         = PixNumOut;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);



    %% Unpack connectivity
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

   
    h = 1e-4;

  
D = zeros(3*N, 3*N);

% Use Parallel Loop (Requires Parallel Computing Toolbox)
for j = 1:3*N

    f_plus = f;
    f_minus = f;

    f_plus(j) = f(j) + h;
    f_minus(j) = f(j) - h;

    S_plus = f_plus(idxS); C_plus = f_plus(idxC); I_plus = f_plus(idxI);
    S_minus = f_minus(idxS); C_minus = f_minus(idxC); I_minus = f_minus(idxI);

    L4E_plus = C_SS_mean*S_plus + C_CS_mean*S_plus + C_IS_mean*S_plus + ...
               C_SC_mean*C_plus + C_CC_mean*C_plus + C_IC_mean*C_plus;

    L4I_plus = C_SI_mean*I_plus + C_CI_mean*I_plus + C_II_mean*I_plus;

    L4E_minus = C_SS_mean*S_minus + C_CS_mean*S_minus + C_IS_mean*S_minus + ...
                C_SC_mean*C_minus + C_CC_mean*C_minus + C_IC_mean*C_minus;

    L4I_minus = C_SI_mean*I_minus + C_CI_mean*I_minus + C_II_mean*I_minus;

    % S Component
    Psi_S_p = local_Psi('S', L4E_plus,  L4I_plus,  PixInptCtgrUse);
    Psi_S_m = local_Psi('S', L4E_minus, L4I_minus, PixInptCtgrUse);
    Diff_S = (Psi_S_p - Psi_S_m) / (2*h);

    % C Component
    Psi_C_p = local_Psi('C', L4E_plus,  L4I_plus,  PixInptCtgrUse);
    Psi_C_m = local_Psi('C', L4E_minus, L4I_minus, PixInptCtgrUse);
    Diff_C = (Psi_C_p - Psi_C_m) / (2*h);

    % I Component
    Psi_I_p = local_Psi('I', L4E_plus,  L4I_plus,  PixInptCtgrUse);
    Psi_I_m = local_Psi('I', L4E_minus, L4I_minus, PixInptCtgrUse);
    Diff_I = (Psi_I_p - Psi_I_m) / (2*h);

    D(:, j) = [Diff_S; Diff_C; Diff_I];

end



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
