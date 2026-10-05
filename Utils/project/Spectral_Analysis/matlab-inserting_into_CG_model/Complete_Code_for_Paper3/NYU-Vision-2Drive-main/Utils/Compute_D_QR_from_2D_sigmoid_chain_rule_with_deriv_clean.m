function J = Compute_D_QR_from_2D_sigmoid_chain_rule_with_deriv_clean(f, ParaODE)
% Analytic Jacobian using *_with_deriv (dPsi/dL4EUse, dPsi/dL4IUse) and explicit
% chain rule: D_QR = diag(dPsi_Q/dL4EUse)*dL4EUse/dR + diag(dPsi_Q/dL4IUse)*dL4IUse/dR.

    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;

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

    % Forward drives (as currently defined in the codebase)
    L4EUse = C_SS_mean*S + C_CS_mean*S + C_IS_mean*S + ...
             C_SC_mean*C + C_CC_mean*C + C_IC_mean*C;

    L4IUse = C_SI_mean*I + C_CI_mean*I + C_II_mean*I;

    L4EUse = L4EUse(:);
    L4IUse = L4IUse(:);

    % Jacobian factors: dL4EUse/dR, dL4IUse/dR
    dL4E_dS = C_SS_mean + C_CS_mean + C_IS_mean;
    dL4E_dC = C_SC_mean + C_CC_mean + C_IC_mean;
    dL4E_dI = 0 * C_II_mean; % explicitly zero matrix of matching size

    dL4I_dS = 0 * C_SI_mean;
    dL4I_dC = 0 * C_CI_mean;
    dL4I_dI = C_SI_mean + C_CI_mean + C_II_mean;

    % Per-pixel derivatives of Psi w.r.t. L4EUse/L4IUse
    [~, dPsiS_dL4E, dPsiS_dL4I] = local_Psi_with_deriv('S', L4EUse, L4IUse, PixInptCtgrUse);
    [~, dPsiC_dL4E, dPsiC_dL4I] = local_Psi_with_deriv('C', L4EUse, L4IUse, PixInptCtgrUse);
    [~, dPsiI_dL4E, dPsiI_dL4I] = local_Psi_with_deriv('I', L4EUse, L4IUse, PixInptCtgrUse);

    % Diagonalize the per-pixel slopes
    diagS_L4E = diag(dPsiS_dL4E);
    diagS_L4I = diag(dPsiS_dL4I);

    diagC_L4E = diag(dPsiC_dL4E);
    diagC_L4I = diag(dPsiC_dL4I);

    diagI_L4E = diag(dPsiI_dL4E);
    diagI_L4I = diag(dPsiI_dL4I);

    % Assemble blocks via chain rule
    D_SS = diagS_L4E * dL4E_dS + diagS_L4I * dL4I_dS;
    D_SC = diagS_L4E * dL4E_dC + diagS_L4I * dL4I_dC;
    D_SI = diagS_L4E * dL4E_dI + diagS_L4I * dL4I_dI;

    D_CS = diagC_L4E * dL4E_dS + diagC_L4I * dL4I_dS;
    D_CC = diagC_L4E * dL4E_dC + diagC_L4I * dL4I_dC;
    D_CI = diagC_L4E * dL4E_dI + diagC_L4I * dL4I_dI;

    D_IS = diagI_L4E * dL4E_dS + diagI_L4I * dL4I_dS;
    D_IC = diagI_L4E * dL4E_dC + diagI_L4I * dL4I_dC;
    D_II = diagI_L4E * dL4E_dI + diagI_L4I * dL4I_dI;

    J = [D_SS, D_SC, D_SI;
         D_CS, D_CC, D_CI;
         D_IS, D_IC, D_II];
end
