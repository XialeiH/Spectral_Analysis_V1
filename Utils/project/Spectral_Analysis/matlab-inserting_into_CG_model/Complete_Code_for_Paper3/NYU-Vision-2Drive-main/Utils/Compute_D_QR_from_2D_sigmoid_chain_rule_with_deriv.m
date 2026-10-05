function J = Compute_D_QR_from_2D_sigmoid_chain_rule_with_deriv(f, ParaODE)
%COMPUTE_D_QR_FROM_2D_SIGMOID_CHAIN_RULE_WITH_DERIV
%   Matches Compute_D_QR_from_2D_sigmoid_chain_rule but obtains the
%   derivatives analytically using local_Psi_with_deriv and the *_with_deriv
%   sigmoid functions, all evaluated at the provided state f.

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

    L4SE = C_SS_mean*S + C_SC_mean*C;
    L4CE = C_CS_mean*S + C_CC_mean*C;
    L4IE = C_IS_mean*S + C_IC_mean*C;

    L4SI = C_SI_mean*I;
    L4CI = C_CI_mean*I;
    L4II = C_II_mean*I;

    L4SEp = mean(L4SE./L4EUse);L4CEp = mean(L4CE./L4EUse);L4IEp = mean(L4IE./L4EUse);
    L4SIp = mean(L4SI./L4IUse);L4CIp = mean(L4CI./L4IUse);L4IIp = mean(L4II./L4IUse);
    


    % evaluate Psi and its derivatives for each cell type (w.r.t. scaled products)
    [~, dPsiS_dL4E, dPsiS_dL4I] = local_Psi_with_deriv('S', L4EUse, L4IUse, PixInptCtgrUse);
    [~, dPsiC_dL4E, dPsiC_dL4I] = local_Psi_with_deriv('C', L4EUse, L4IUse, PixInptCtgrUse);
    [~, dPsiI_dL4E, dPsiI_dL4I] = local_Psi_with_deriv('I', L4EUse, L4IUse, PixInptCtgrUse);

    % build diagonal blocks directly from analytic derivatives
    diagS_L4E = diag(dPsiS_dL4E);
    diagS_L4I = diag(dPsiS_dL4I);

    diagC_L4E = diag(dPsiC_dL4E);
    diagC_L4I = diag(dPsiC_dL4I);

    diagI_L4E = diag(dPsiI_dL4E);
    diagI_L4I = diag(dPsiI_dL4I);

    D_SS = diagS_L4E * C_SS_mean;
    D_SC = diagS_L4E * C_SC_mean;
    D_SI = diagS_L4I * C_SI_mean;

    D_CS = diagC_L4E * C_CS_mean;
    D_CC = diagC_L4E * C_CC_mean;
    D_CI = diagC_L4I * C_CI_mean;

    D_IS = diagI_L4E * C_IS_mean;
    D_IC = diagI_L4E * C_IC_mean;
    D_II = diagI_L4I * C_II_mean;


    J = [D_SS,D_SC,D_SI; D_CS,D_CC,D_CI;D_IS,D_IC,D_II];
end
