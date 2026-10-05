

function [D_SS,D_SC,D_SI, ...
          D_CS,D_CC,D_CI, ...
          D_IS,D_IC,D_II] = ...
    Compute_D_QR_from_2D_sigmoid(f, ParaODE)


    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
    N         = PixNumOut;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

    S = f(idxS);
    C = f(idxC);
    I = f(idxI);

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

    %% Build the (shared) drives L4EUse and L4IUse at this operating point

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

    %% Finite-difference step
    h = 1e-3;

    %% ---------------- Q = S ----------------

    % ----- D_SS (Q=S, R=S → perturb L4E, C_SS_mean) -----
    D_SS = zeros(N, N);
    for j = 1:N
        for i = 1:N
        L4E_plus  = L4EUse;
        L4E_minus = L4EUse;
        L4I_plus  = L4IUse;
        L4I_minus = L4IUse;

        % R = S: perturb excitatory drive at pixel i
        L4E_plus(j)  = L4E_plus(j)  + h;
        L4E_minus(j) = L4E_minus(j) - h;

        Psi_Q_plus  = local_Psi('S', L4E_plus,  L4I_plus,  PixInptCtgrUse);
        Psi_Q_minus = local_Psi('S', L4E_minus, L4I_minus, PixInptCtgrUse);

        slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

        D_SS(i,j) = slope_ij * C_SS_mean(i,j);
        
        end
    end

    % ----- D_SC (Q=S, R=C, use C_SC_mean) -----
    D_SC = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = C: still excitatory drive at pixel j
            L4E_plus(j)  = L4E_plus(j)  + h;
            L4E_minus(j) = L4E_minus(j) - h;

            Psi_Q_plus  = local_Psi('S', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('S', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_SC(i,j) = slope_ij * C_SC_mean(i,j);
        end
    end

    % ----- D_SI (Q=S, R=I, use C_SI_mean) -----
    D_SI = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = I: perturb inhibitory drive at pixel j
            L4I_plus(j)  = L4I_plus(j)  + h;
            L4I_minus(j) = L4I_minus(j) - h;

            Psi_Q_plus  = local_Psi('S', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('S', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_SI(i,j) = slope_ij * C_SI_mean(i,j);
        end
    end

    %% ---------------- Q = C ----------------

    % ----- D_CS (Q=C, R=S, use C_CS_mean) -----
    D_CS = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = S: perturb excitatory drive at pixel j
            L4E_plus(j)  = L4E_plus(j)  + h;
            L4E_minus(j) = L4E_minus(j) - h;

            Psi_Q_plus  = local_Psi('C', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('C', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_CS(i,j) = slope_ij * C_CS_mean(i,j);
        end
    end

    % ----- D_CC (Q=C, R=C, use C_CC_mean) -----
    D_CC = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = C: excitatory drive at pixel j
            L4E_plus(j)  = L4E_plus(j)  + h;
            L4E_minus(j) = L4E_minus(j) - h;

            Psi_Q_plus  = local_Psi('C', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('C', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_CC(i,j) = slope_ij * C_CC_mean(i,j);
        end
    end

    % ----- D_CI (Q=C, R=I, use C_CI_mean) -----
    D_CI = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = I: perturb inhibitory drive at pixel j
            L4I_plus(j)  = L4I_plus(j)  + h;
            L4I_minus(j) = L4I_minus(j) - h;

            Psi_Q_plus  = local_Psi('C', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('C', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_CI(i,j) = slope_ij * C_CI_mean(i,j);
        end
    end

    %% ---------------- Q = I ----------------

    % ----- D_IS (Q=I, R=S, use C_IS_mean) -----
    D_IS = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = S: perturb excitatory drive at pixel j
            L4E_plus(j)  = L4E_plus(j)  + h;
            L4E_minus(j) = L4E_minus(j) - h;

            Psi_Q_plus  = local_Psi('I', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('I', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_IS(i,j) = slope_ij * C_IS_mean(i,j);
        end
    end

    % ----- D_IC (Q=I, R=C, use C_IC_mean) -----
    D_IC = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = C: excitatory drive at pixel j
            L4E_plus(j)  = L4E_plus(j)  + h;
            L4E_minus(j) = L4E_minus(j) - h;

            Psi_Q_plus  = local_Psi('I', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('I', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_IC(i,j) = slope_ij * C_IC_mean(i,j);
        end
    end

    % ----- D_II (Q=I, R=I, use C_II_mean) -----
    D_II = zeros(N, N);
    for j = 1:N
        for i = 1:N
            L4E_plus  = L4EUse;
            L4E_minus = L4EUse;
            L4I_plus  = L4IUse;
            L4I_minus = L4IUse;

            % R = I: perturb inhibitory drive at pixel j
            L4I_plus(j)  = L4I_plus(j)  + h;
            L4I_minus(j) = L4I_minus(j) - h;

            Psi_Q_plus  = local_Psi('I', L4E_plus,  L4I_plus,  PixInptCtgrUse);
            Psi_Q_minus = local_Psi('I', L4E_minus, L4I_minus, PixInptCtgrUse);

            slope_ij = (Psi_Q_plus(i) - Psi_Q_minus(i)) / (2*h);

            D_II(i,j) = slope_ij * C_II_mean(i,j);
        end
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

% 
% function [D_SS,D_SC,D_SI, ...
%           D_CS,D_CC,D_CI, ...
%           D_IS,D_IC,D_II] = ...
%     Compute_D_QR_from_2D_sigmoid(f, ParaODE)
% % Compute_D_QR_from_2D_sigmoid
% %
% % For each postsynaptic population Q ∈ {S,C,I} and each presynaptic
% % population R ∈ {S,C,I}, compute the full matrix
% %
% %    D_QR(i,j) = ∂Ψ_Q(i) / ∂(R-input at pixel j)
% %
% % where Ψ_Q is the local response function built from the 16 LGN–L6
% % libraries and the pixelwise PixInptCtgrUse (size N×4×4).
% %
% % This code:
% %   1) builds the *drives* L4EUse and L4IUse from the operating point f
% %      using ParaODE.C_*_mean,
% %   2) for each Q, uses finite differences in L4EUse and L4IUse at each
% %      pixel to obtain the local slopes dΨ_Q(i)/dE(i), dΨ_Q(i)/dI(i),
% %      with the pixel-specific PixInptCtgrUse included,
% %   3) uses the chain rule with the connectivity matrices C_*_mean to
% %      form full N×N matrices D_QR.
% %
% % IMPORTANT:
% %   - D_QR here are already full matrices. If elsewhere you previously
% %     did   J11 = D_SS * C_SS_mean, etc., you should now use
% %     J11 = D_SS, J12 = D_SC, etc., because the connectivity has been
% %     folded into D_QR.
% 
%     %% Basic sizes and indexing
%     PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
% 
%     idxS = 1:PixNumOut;
%     idxC = PixNumOut + (1:PixNumOut);
%     idxI = 2*PixNumOut + (1:PixNumOut);
% 
%     S = f(idxS);
%     C = f(idxC);
%     I = f(idxI);
% 
%     %% Unpack connectivity and category weights
%     C_SS_mean = ParaODE.C_SS_mean;
%     C_CS_mean = ParaODE.C_CS_mean;
%     C_IS_mean = ParaODE.C_IS_mean;
% 
%     C_SC_mean = ParaODE.C_SC_mean;
%     C_CC_mean = ParaODE.C_CC_mean;
%     C_IC_mean = ParaODE.C_IC_mean;
% 
%     C_SI_mean = ParaODE.C_SI_mean;
%     C_CI_mean = ParaODE.C_CI_mean;
%     C_II_mean = ParaODE.C_II_mean;
% 
%     PixInptCtgrUse = ParaODE.PixInptCtgrUse;  % N×4×4
% 
%     %% Build the (shared) drives L4EUse and L4IUse at this operating point
%     % NOTE: This follows your existing formula; each pixel has one
%     % excitatory drive L4EUse(i) and one inhibitory drive L4IUse(i).
% 
%     L4EUse = C_SS_mean*S + ...
%              C_CS_mean*S + ...
%              C_SI_mean*I + ...
%              C_SC_mean*C + ...
%              C_CC_mean*C + ...
%              C_CI_mean*I;
% 
%     L4IUse = C_SI_mean*I + ...
%              C_CI_mean*I + ...
%              C_II_mean*I;
% 
%     L4EUse = L4EUse(:);
%     L4IUse = L4IUse(:);
% 
%     %% Choose finite-difference steps based on the scale of drives
%     scaleE = max(1, median(abs(L4EUse)));
%     scaleI = max(1, median(abs(L4IUse)));
%     hE = 1e-2 * scaleE;     % ~1% of typical excitatory drive
%     hI = 1e-2 * scaleI;     % ~1% of typical inhibitory drive
% 
%     %% 1) Local derivatives for each population Q, using the *inlined* library
% 
%     [dS_dE, dS_dI] = local_derivatives_for_Q('S', ...
%                                              L4EUse, L4IUse, ...
%                                              PixInptCtgrUse, ...
%                                              hE, hI, PixNumOut);
% 
%     [dC_dE, dC_dI] = local_derivatives_for_Q('C', ...
%                                              L4EUse, L4IUse, ...
%                                              PixInptCtgrUse, ...
%                                              hE, hI, PixNumOut);
% 
%     [dI_dE, dI_dI] = local_derivatives_for_Q('I', ...
%                                              L4EUse, L4IUse, ...
%                                              PixInptCtgrUse, ...
%                                              hE, hI, PixNumOut);
% 
%     % Diagonal local gains wrt E and I
%     G_S_E = spdiags(dS_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_S/∂E (local)
%     G_S_I = spdiags(dS_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_S/∂I (local)
% 
%     G_C_E = spdiags(dC_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_C/∂E
%     G_C_I = spdiags(dC_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_C/∂I
% 
%     G_I_E = spdiags(dI_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_I/∂E
%     G_I_I = spdiags(dI_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_I/∂I
% 
%     %% 2) Chain rule: fold connectivity into D_QR so that
%     %      D_QR(i,j) = ∂Ψ_Q(i)/∂(R-input at pixel j)
%     %
%     % Here we follow your Jacobian-block pattern:
%     %   originally: J11 = G_S_E*C_SS_mean; J12 = G_S_E*C_CS_mean; ...
%     % so now we define D_QR directly as these full products.
% 
%     % Postsynaptic S:
%     D_SS = G_S_E * C_SS_mean;   % ∂Ψ_S / ∂S-input
%     D_SC = G_S_E * C_CS_mean;   % ∂Ψ_S / ∂C-input
%     D_SI = G_S_I * C_IS_mean;   % ∂Ψ_S / ∂I-input
% 
%     % Postsynaptic C:
%     D_CS = G_C_E * C_SC_mean;   % ∂Ψ_C / ∂S-input
%     D_CC = G_C_E * C_CC_mean;   % ∂Ψ_C / ∂C-input
%     D_CI = G_C_I * C_IC_mean;   % ∂Ψ_C / ∂I-input
% 
%     % Postsynaptic I:
%     D_IS = G_I_E * C_SI_mean;   % ∂Ψ_I / ∂S-input
%     D_IC = G_I_E * C_CI_mean;   % ∂Ψ_I / ∂C-input
%     D_II = G_I_I * C_II_mean;   % ∂Ψ_I / ∂I-input
% end
% 
% %======================================================================
% % Local helper: compute dΨ_Q/dE and dΨ_Q/dI per pixel for a given Q,
% % using finite differences on the drives and the inlined 16-library map.
% %======================================================================
% function [dPsi_dE, dPsi_dI] = local_derivatives_for_Q( ...
%             celltype, L4EUse, L4IUse, PixInptCtgrUse, hE, hI, PixNumOut)
% 
%     dPsi_dE = zeros(PixNumOut,1);
%     dPsi_dI = zeros(PixNumOut,1);
% 
%     for i = 1:PixNumOut
%         e_i    = zeros(PixNumOut,1);
%         e_i(i) = 1;
% 
%         % --- perturb excitatory drive at pixel i ---
%         E_plus  = L4EUse + hE * e_i;
%         E_minus = L4EUse - hE * e_i;
% 
%         Psi_plus_E  = local_Psi(celltype, E_plus,  L4IUse, PixInptCtgrUse);
%         Psi_minus_E = local_Psi(celltype, E_minus, L4IUse, PixInptCtgrUse);
% 
%         dPsi_dE(i) = (Psi_plus_E(i) - Psi_minus_E(i)) / (2*hE);
% 
%         % --- perturb inhibitory drive at pixel i ---
%         I_plus  = L4IUse + hI * e_i;
%         I_minus = L4IUse - hI * e_i;
% 
%         Psi_plus_I  = local_Psi(celltype, L4EUse, I_plus,  PixInptCtgrUse);
%         Psi_minus_I = local_Psi(celltype, L4EUse, I_minus, PixInptCtgrUse);
% 
%         dPsi_dI(i) = (Psi_plus_I(i) - Psi_minus_I(i)) / (2*hI);
%     end
% end
% 
% %======================================================================
% % Inlined local response function:
% %   Ψ_Q = weighted sum over 16 LGN–L6 library outputs, pixelwise.
% %
% % This replaces LDEIterFunc_Grating_16Func_noFor_0324, but uses the same
% % "essence": 16 library functions + pixelwise weights PixInptCtgrUse.
% %======================================================================








% 
% function [D_SS,D_SC,D_SI, ...
%           D_CS,D_CC,D_CI, ...
%           D_IS,D_IC,D_II] = ...
%     Compute_D_QR_from_2D_sigmoid(f, ParaODE)
% % Compute_D_QR_from_2D_sigmoid
% %
% % 1. Aggregates inputs from S and C into a single EXCITATORY Drive (TotalE).
% % 2. Aggregates inputs from I into a single INHIBITORY Drive (TotalI).
% % 3. Computes the LOCAL GAINS (Slope of Sigmoid) by perturbing TotalE and TotalI.
% %    (Satisfies the requirement that Sigmoids need full Vector inputs).
% % 4. Applies Chain Rule:
% %    - Change in S or C -> affects TotalE -> multiplied by G_E
% %    - Change in I      -> affects TotalI -> multiplied by G_I
% 
%     %% Basic sizes and indexing
%     PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY; % e.g., 1600
%     N = PixNumOut;
%     idxS = 1:N;
%     idxC = N + (1:N);
%     idxI = 2*N + (1:N);
% 
%     %% Unpack Connectivity
%     C_SS_mean = ParaODE.C_SS_mean;
%     C_CS_mean = ParaODE.C_CS_mean;
%     C_IS_mean = ParaODE.C_IS_mean;
% 
%     C_SC_mean = ParaODE.C_SC_mean;
%     C_CC_mean = ParaODE.C_CC_mean;
%     C_IC_mean = ParaODE.C_IC_mean;
% 
%     C_SI_mean = ParaODE.C_SI_mean;
%     C_CI_mean = ParaODE.C_CI_mean;
%     C_II_mean = ParaODE.C_II_mean;
% 
%     PixInptCtgrUse = ParaODE.PixInptCtgrUse;
% 
%     %% 1) Calculate Total Drives (The "Operating Point")
%     % IMPORTANT: S and C are put together here to form the Excitatory Drive.
% 
%     S = f(idxS);
%     C = f(idxC);
%     I = f(idxI);
% 
%     % --- For Postsynaptic S ---
%     % S receives Excitatory input from both S and C
%     TotalE_S = C_SS_mean*S + C_SC_mean*C;   
%     % S receives Inhibitory input from I
%     TotalI_S = C_SI_mean*I;                 
% 
%     % --- For Postsynaptic C ---
%     % C receives Excitatory input from both S and C
%     TotalE_C = C_CS_mean*S + C_CC_mean*C;   
%     TotalI_C = C_CI_mean*I;                 
% 
%     % --- For Postsynaptic I ---
%     % I receives Excitatory input from both S and C
%     TotalE_I = C_IS_mean*S + C_IC_mean*C;   
%     TotalI_I = C_II_mean*I;                 
% 
%     %% 2) Choose Finite-Difference Steps (h) based on Drive Magnitude
%     scaleE = max(1, median(abs([TotalE_S; TotalE_C; TotalE_I]))); 
%     scaleI = max(1, median(abs([TotalI_S; TotalI_C; TotalI_I])));
%     hE = 1e-4 * scaleE;
%     hI = 1e-4 * scaleI;
% 
%     %% 3) Compute Local Gains (Vectorized)
%     % We compute the slope of the sigmoid w.r.t the Total Excitatory Drive
%     % and the Total Inhibitory Drive.
% 
%     % Gain for S (Slope of Psi_S)
%     [dS_dTotalE, dS_dTotalI] = get_vectorized_gains('S', TotalE_S, TotalI_S, PixInptCtgrUse, hE, hI);
% 
%     % Gain for C (Slope of Psi_C)
%     [dC_dTotalE, dC_dTotalI] = get_vectorized_gains('C', TotalE_C, TotalI_C, PixInptCtgrUse, hE, hI);
% 
%     % Gain for I (Slope of Psi_I)
%     [dI_dTotalE, dI_dTotalI] = get_vectorized_gains('I', TotalE_I, TotalI_I, PixInptCtgrUse, hE, hI);
% 
%     %% 4) Create Diagonal Gain Matrices
%     % G_X_E answers: "How much does FiringRate X change if Total E-input changes?"
% 
%     G_S_E = spdiags(dS_dTotalE, 0, N, N);  
%     G_S_I = spdiags(dS_dTotalI, 0, N, N);
% 
%     G_C_E = spdiags(dC_dTotalE, 0, N, N); 
%     G_C_I = spdiags(dC_dTotalI, 0, N, N);
% 
%     G_I_E = spdiags(dI_dTotalE, 0, N, N); 
%     G_I_I = spdiags(dI_dTotalI, 0, N, N);
% 
%     %% 5) Chain Rule: D_QR = Gain * Connectivity
%     % Note: G_S_E applies to BOTH D_SS and D_SC, because S and C both 
%     % contribute to the Excitatory drive.
% 
%     % Postsynaptic S
%     D_SS = G_S_E * C_SS_mean;   % Input from S -> E-Drive -> Output S
%     D_SC = G_S_E * C_SC_mean;   % Input from C -> E-Drive -> Output S
%     D_SI = G_S_I * C_SI_mean;   % Input from I -> I-Drive -> Output S
% 
%     % Postsynaptic C
%     D_CS = G_C_E * C_CS_mean;
%     D_CC = G_C_E * C_CC_mean;
%     D_CI = G_C_I * C_CI_mean;
% 
%     % Postsynaptic I
%     D_IS = G_I_E * C_IS_mean;
%     D_IC = G_I_E * C_IC_mean;
%     D_II = G_I_I * C_II_mean;
% end
% 
% %======================================================================
% % Helper: Compute Gains for ALL pixels simultaneously
% % Returns the derivative w.r.t TotalE and TotalI
% %======================================================================
% function [dPsi_dE, dPsi_dI] = get_vectorized_gains(celltype, E_vec, I_vec, PixInptCtgrUse, hE, hI)
% 
%     % --- Perturb Total Excitatory Drive (Global Shift) ---
%     % This captures sensitivity to the COMBINED S+C input
%     Psi_plus_E  = local_Psi_vector(celltype, E_vec + hE, I_vec,      PixInptCtgrUse);
%     Psi_minus_E = local_Psi_vector(celltype, E_vec - hE, I_vec,      PixInptCtgrUse);
%     dPsi_dE     = (Psi_plus_E - Psi_minus_E) / (2*hE);
% 
%     % --- Perturb Total Inhibitory Drive (Global Shift) ---
%     Psi_plus_I  = local_Psi_vector(celltype, E_vec,      I_vec + hI, PixInptCtgrUse);
%     Psi_minus_I = local_Psi_vector(celltype, E_vec,      I_vec - hI, PixInptCtgrUse);
%     dPsi_dI     = (Psi_plus_I - Psi_minus_I) / (2*hI);
% end
% 
% %======================================================================
% % Inlined local response function: 16-Library Sigmoid Map
% % Accepts FULL VECTORS (Nx1) to satisfy reshape() requirements
% %======================================================================
% function Psi_Q = local_Psi_vector(celltype, L4EUse, L4IUse, PixInptCtgrUse)
% 
%     Npix = numel(L4EUse);
%     LibyAll = zeros(Npix,4,4);  % [pixel, LGN, L6]
% 
%     switch celltype
%         case 'S'
%             % LGN c1
%             LibyAll(:,1,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = S_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = S_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = S_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c2
%             LibyAll(:,2,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = S_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = S_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = S_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c3
%             LibyAll(:,3,1) = S_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = S_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = S_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = S_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c4
%             LibyAll(:,4,1) = S_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = S_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = S_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = S_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
% 
%         case 'C'
%             LibyAll(:,1,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = C_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = C_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = C_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = C_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = C_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = C_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,1) = C_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = C_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = C_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = C_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,1) = C_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = C_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = C_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = C_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
% 
%         case 'I'
%             LibyAll(:,1,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = I_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = I_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = I_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = I_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = I_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = I_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,1) = I_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = I_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = I_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = I_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,1) = I_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = I_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = I_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = I_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
%     end
% 
%     Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2,3]);
%     Psi_Q = Psi_Q(:);  
% end









% function [D_SS,D_SC,D_SI, ...
%           D_CS,D_CC,D_CI, ...
%           D_IS,D_IC,D_II] = ...
%     Compute_D_QR_from_2D_sigmoid(f, ParaODE)
% % Compute_D_QR_from_2D_sigmoid
% %
% %  Explicitly compute the 9 Jacobian blocks D_PQ by finite differences.
% %  P,Q ∈ {S,C,I}, each block is N×N.
% %
% %  For each presynaptic population Q and each neuron j:
% %     • Perturb f_Q(j) by ±h_Q
% %     • Recompute all 3 postsynaptic outputs Ψ_S, Ψ_C, Ψ_I
% %     • Use central difference to fill the j-th column of
% %           D_SQ(:,j), D_CQ(:,j), D_IQ(:,j)
% %
% %  This treats the local response functions Ψ_Q as depending on the
% %  FULL vectors of excitatory and inhibitory drives, so the resulting
% %  D_PQ are fully populated (not diagonal).
% 
%     %% Basic sizes and indexing
%     PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY; % e.g., 1600
%     N = PixNumOut;
%     idxS = 1:N;
%     idxC = N + (1:N);
%     idxI = 2*N + (1:N);
% 
%     %% Unpack f into S, C, I
%     S = f(idxS);
%     C = f(idxC);
%     I = f(idxI);
% 
%     %% Unpack Connectivity
%     C_SS_mean = ParaODE.C_SS_mean;
%     C_CS_mean = ParaODE.C_CS_mean;
%     C_IS_mean = ParaODE.C_IS_mean;
% 
%     C_SC_mean = ParaODE.C_SC_mean;
%     C_CC_mean = ParaODE.C_CC_mean;
%     C_IC_mean = ParaODE.C_IC_mean;
% 
%     C_SI_mean = ParaODE.C_SI_mean;
%     C_CI_mean = ParaODE.C_CI_mean;
%     C_II_mean = ParaODE.C_II_mean;
% 
%     PixInptCtgrUse = ParaODE.PixInptCtgrUse;
% 
%     %% Choose finite-difference steps for S, C, I
%     scaleS = max(1, median(abs(S)));
%     scaleC = max(1, median(abs(C)));
%     scaleI = max(1, median(abs(I)));
% 
%     hS = 1e-4 * scaleS;
%     hC = 1e-4 * scaleC;
%     hI = 1e-4 * scaleI;
% 
%     %% Allocate the 9 N×N blocks
%     D_SS = zeros(N,N);
%     D_SC = zeros(N,N);
%     D_SI = zeros(N,N);
% 
%     D_CS = zeros(N,N);
%     D_CC = zeros(N,N);
%     D_CI = zeros(N,N);
% 
%     D_IS = zeros(N,N);
%     D_IC = zeros(N,N);
%     D_II = zeros(N,N);
% 
%     %% ===== 1) Presynaptic S: fill D_SS, D_CS, D_IS =====
%     for j = 1:N
%         % Perturb S(j) up and down
%         S_plus  = S; S_plus(j)  = S_plus(j)  + hS;
%         S_minus = S; S_minus(j) = S_minus(j) - hS;
% 
%         [Psi_S_plus, Psi_C_plus, Psi_I_plus] = ...
%             compute_Psi_all(S_plus, C, I, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         [Psi_S_minus, Psi_C_minus, Psi_I_minus] = ...
%             compute_Psi_all(S_minus, C, I, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         D_SS(:,j) = (Psi_S_plus - Psi_S_minus) / (2*hS);  % ∂Ψ_S / ∂S_j
%         D_CS(:,j) = (Psi_C_plus - Psi_C_minus) / (2*hS);  % ∂Ψ_C / ∂S_j
%         D_IS(:,j) = (Psi_I_plus - Psi_I_minus) / (2*hS);  % ∂Ψ_I / ∂S_j
%     end
% 
%     %% ===== 2) Presynaptic C: fill D_SC, D_CC, D_IC =====
%     for j = 1:N
%         % Perturb C(j) up and down
%         C_plus  = C; C_plus(j)  = C_plus(j)  + hC;
%         C_minus = C; C_minus(j) = C_minus(j) - hC;
% 
%         [Psi_S_plus, Psi_C_plus, Psi_I_plus] = ...
%             compute_Psi_all(S, C_plus, I, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         [Psi_S_minus, Psi_C_minus, Psi_I_minus] = ...
%             compute_Psi_all(S, C_minus, I, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         D_SC(:,j) = (Psi_S_plus - Psi_S_minus) / (2*hC);  % ∂Ψ_S / ∂C_j
%         D_CC(:,j) = (Psi_C_plus - Psi_C_minus) / (2*hC);  % ∂Ψ_C / ∂C_j
%         D_IC(:,j) = (Psi_I_plus - Psi_I_minus) / (2*hC);  % ∂Ψ_I / ∂C_j
%     end
% 
%     %% ===== 3) Presynaptic I: fill D_SI, D_CI, D_II =====
%     for j = 1:N
%         % Perturb I(j) up and down
%         I_plus  = I; I_plus(j)  = I_plus(j)  + hI;
%         I_minus = I; I_minus(j) = I_minus(j) - hI;
% 
%         [Psi_S_plus, Psi_C_plus, Psi_I_plus] = ...
%             compute_Psi_all(S, C, I_plus, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         [Psi_S_minus, Psi_C_minus, Psi_I_minus] = ...
%             compute_Psi_all(S, C, I_minus, ...
%                 C_SS_mean,C_SC_mean,C_SI_mean, ...
%                 C_CS_mean,C_CC_mean,C_CI_mean, ...
%                 C_IS_mean,C_IC_mean,C_II_mean, ...
%                 PixInptCtgrUse);
% 
%         D_SI(:,j) = (Psi_S_plus - Psi_S_minus) / (2*hI);  % ∂Ψ_S / ∂I_j
%         D_CI(:,j) = (Psi_C_plus - Psi_C_minus) / (2*hI);  % ∂Ψ_C / ∂I_j
%         D_II(:,j) = (Psi_I_plus - Psi_I_minus) / (2*hI);  % ∂Ψ_I / ∂I_j
%     end
% end
% 
% %======================================================================
% % Helper: given S, C, I, compute all three outputs Ψ_S, Ψ_C, Ψ_I
% %======================================================================
% function [Psi_S, Psi_C, Psi_I] = compute_Psi_all( ...
%         S, C, I, ...
%         C_SS_mean,C_SC_mean,C_SI_mean, ...
%         C_CS_mean,C_CC_mean,C_CI_mean, ...
%         C_IS_mean,C_IC_mean,C_II_mean, ...
%         PixInptCtgrUse)
% 
%     % --- Postsynaptic S drives ---
%     L4EUse = C_SS_mean*S + ...
%              C_CS_mean*S + ...
%              C_SI_mean*I + ...
%              C_SC_mean*C + ...
%              C_CC_mean*C + ...
%              C_CI_mean*I;
% 
%     L4IUse = C_SI_mean*I + ...
%              C_CI_mean*I + ...
%              C_II_mean*I;
% 
%     L4EUse = L4EUse(:);
%     L4IUse = L4IUse(:);
% 
%     % Local nonlinear maps (use your 16-function library)
%     Psi_S = local_Psi_vector('S', L4EUse, L4IUse, PixInptCtgrUse);
%     Psi_C = local_Psi_vector('C', L4EUse, L4IUse, PixInptCtgrUse);
%     Psi_I = local_Psi_vector('I', L4EUse, L4IUse, PixInptCtgrUse);
% end
% 
% %======================================================================
% % Inlined local response function: 16-Library Sigmoid Map
% % Accepts FULL VECTORS (Nx1) to satisfy reshape() requirements
% %======================================================================
% function Psi_Q = local_Psi_vector(celltype, L4EUse, L4IUse, PixInptCtgrUse)
% 
%     Npix = numel(L4EUse);
%     LibyAll = zeros(Npix,4,4);  % [pixel, LGN, L6]
% 
%     switch celltype
%         case 'S'
%             % LGN c1
%             LibyAll(:,1,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = S_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = S_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = S_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c2
%             LibyAll(:,2,1) = S_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = S_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = S_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = S_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c3
%             LibyAll(:,3,1) = S_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = S_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = S_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = S_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             % LGN c4
%             LibyAll(:,4,1) = S_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = S_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = S_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = S_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
% 
%         case 'C'
%             LibyAll(:,1,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = C_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = C_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = C_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,1) = C_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = C_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = C_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = C_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,1) = C_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = C_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = C_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = C_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,1) = C_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = C_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = C_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = C_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
% 
%         case 'I'
%             LibyAll(:,1,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,2) = I_LGNc1_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,3) = I_LGNc1_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,1,4) = I_LGNc1_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,1) = I_LGNc1_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,2) = I_LGNc2_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,3) = I_LGNc2_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,2,4) = I_LGNc2_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,1) = I_LGNc3_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,2) = I_LGNc3_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,3) = I_LGNc3_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,3,4) = I_LGNc3_L6c4_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,1) = I_LGNc4_L6c1_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,2) = I_LGNc4_L6c2_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,3) = I_LGNc4_L6c3_Sigmoid(L4EUse, L4IUse);
%             LibyAll(:,4,4) = I_LGNc4_L6c4_Sigmoid(L4EUse, L4IUse);
%     end
% 
%     Psi_Q = sum(PixInptCtgrUse .* LibyAll, [2,3]);
%     Psi_Q = Psi_Q(:);  
% end