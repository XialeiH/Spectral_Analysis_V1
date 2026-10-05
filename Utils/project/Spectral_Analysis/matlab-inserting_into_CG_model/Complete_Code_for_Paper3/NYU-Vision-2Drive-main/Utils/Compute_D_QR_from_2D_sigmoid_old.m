
function [D_SS,D_SC,D_SI, ...
          D_CS,D_CC,D_CI, ...
          D_IS,D_IC,D_II] = ...
    Compute_D_QR_from_2D_sigmoid_old(f, ParaODE)
% Compute_D_QR_from_2D_sigmoid
%
% For each postsynaptic population Q ∈ {S,C,I} and each presynaptic
% population R ∈ {S,C,I}, compute the full matrix
%
%    D_QR(i,j) = ∂Ψ_Q(i) / ∂(R-input at pixel j)
%
% where Ψ_Q is the local response function built from the 16 LGN–L6
% libraries and the pixelwise PixInptCtgrUse (size N×4×4).
%
% This code:
%   1) builds the *drives* L4EUse and L4IUse from the operating point f
%      using ParaODE.C_*_mean,
%   2) for each Q, uses finite differences in L4EUse and L4IUse at each
%      pixel to obtain the local slopes dΨ_Q(i)/dE(i), dΨ_Q(i)/dI(i),
%      with the pixel-specific PixInptCtgrUse included,
%   3) uses the chain rule with the connectivity matrices C_*_mean to
%      form full N×N matrices D_QR.
%
    %% Basic sizes and indexing
    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

    S = f(idxS);
    C = f(idxC);
    I = f(idxI);

    %% Unpack connectivity and category weights
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
    % NOTE: This follows your existing formula; each pixel has one
    % excitatory drive L4EUse(i) and one inhibitory drive L4IUse(i).

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

    %% Choose finite-difference steps based on the scale of drives
    scaleE = max(1, median(abs(L4EUse)));
    scaleI = max(1, median(abs(L4IUse)));
    hE = 1e-2 * scaleE;     % ~1% of typical excitatory drive
    hI = 1e-2 * scaleI;     % ~1% of typical inhibitory drive

    %% 1) Local derivatives for each population Q, using the *inlined* library

    [dS_dE, dS_dI] = local_derivatives_for_Q('S', ...
                                             L4EUse, L4IUse, ...
                                             PixInptCtgrUse, ...
                                             hE, hI, PixNumOut);

    [dC_dE, dC_dI] = local_derivatives_for_Q('C', ...
                                             L4EUse, L4IUse, ...
                                             PixInptCtgrUse, ...
                                             hE, hI, PixNumOut);

    [dI_dE, dI_dI] = local_derivatives_for_Q('I', ...
                                             L4EUse, L4IUse, ...
                                             PixInptCtgrUse, ...
                                             hE, hI, PixNumOut);

    % Diagonal local gains wrt E and I
    G_S_E = spdiags(dS_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_S/∂E (local)
    G_S_I = spdiags(dS_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_S/∂I (local)

    G_C_E = spdiags(dC_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_C/∂E
    G_C_I = spdiags(dC_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_C/∂I

    G_I_E = spdiags(dI_dE, 0, PixNumOut, PixNumOut);  % ∂Ψ_I/∂E
    G_I_I = spdiags(dI_dI, 0, PixNumOut, PixNumOut);  % ∂Ψ_I/∂I

    %% 2) Chain rule: fold connectivity into D_QR so that
    %      D_QR(i,j) = ∂Ψ_Q(i)/∂(R-input at pixel j)
    %
    % Here we follow your Jacobian-block pattern:
    %   originally: J11 = G_S_E*C_SS_mean; J12 = G_S_E*C_CS_mean; ...
    % so now we define D_QR directly as these full products.

    % Postsynaptic S:
    D_SS = G_S_E ;   % ∂Ψ_S / ∂S-input
    D_SC = G_S_E ;   % ∂Ψ_S / ∂C-input
    D_SI = G_S_I ;   % ∂Ψ_S / ∂I-input

    % Postsynaptic C:
    D_CS = G_C_E ;   % ∂Ψ_C / ∂S-input
    D_CC = G_C_E ;   % ∂Ψ_C / ∂C-input
    D_CI = G_C_I ;   % ∂Ψ_C / ∂I-input

    % Postsynaptic I:
    D_IS = G_I_E ;   % ∂Ψ_I / ∂S-input
    D_IC = G_I_E ;   % ∂Ψ_I / ∂C-input
    D_II = G_I_I ;   % ∂Ψ_I / ∂I-input
end

%======================================================================
% Local helper: compute dΨ_Q/dE and dΨ_Q/dI per pixel for a given Q,
% using finite differences on the drives and the inlined 16-library map.
%======================================================================
function [dPsi_dE, dPsi_dI] = local_derivatives_for_Q( ...
            celltype, L4EUse, L4IUse, PixInptCtgrUse, hE, hI, PixNumOut)

    dPsi_dE = zeros(PixNumOut,1);
    dPsi_dI = zeros(PixNumOut,1);

    for i = 1:PixNumOut
        e_i    = zeros(PixNumOut,1);
        e_i(i) = 1;

        % --- perturb excitatory drive at pixel i ---
        E_plus  = L4EUse + hE * e_i;
        E_minus = L4EUse - hE * e_i;

        Psi_plus_E  = local_Psi(celltype, E_plus,  L4IUse, PixInptCtgrUse);
        Psi_minus_E = local_Psi(celltype, E_minus, L4IUse, PixInptCtgrUse);

        dPsi_dE(i) = (Psi_plus_E(i) - Psi_minus_E(i)) / (2*hE);

        % --- perturb inhibitory drive at pixel i ---
        I_plus  = L4IUse + hI * e_i;
        I_minus = L4IUse - hI * e_i;

        Psi_plus_I  = local_Psi(celltype, L4EUse, I_plus,  PixInptCtgrUse);
        Psi_minus_I = local_Psi(celltype, L4EUse, I_minus, PixInptCtgrUse);

        dPsi_dI(i) = (Psi_plus_I(i) - Psi_minus_I(i)) / (2*hI);
    end
end

%======================================================================
% Inlined local response function:
%   Ψ_Q = weighted sum over 16 LGN–L6 library outputs, pixelwise.
%
% This replaces LDEIterFunc_Grating_16Func_noFor_0324, but uses the same
% "essence": 16 library functions + pixelwise weights PixInptCtgrUse.
%======================================================================

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
