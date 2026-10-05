% Do not need LDEIteration_16FuncMain_CombDom_0324.m by Qingyu!
% Use the same idea but directly out put the corresponding f


function Phi_f = Phi_from_2D_sigmoid(f, ParaODE)
% Phi_from_2D_sigmoid
%   Compute Φ(f) using the 3D local response library (16 functions),
%   inlined from LDEIteration_16FuncMain_CombDom_0324 for ONE step,
%   WITHOUT the mixing step (so this is the pure nonlinear map Φ).
%
%   f is a 3*PixNumOut vector: [S; C; I]
%   ParaODE must contain:
%     N_HCOut, NPixX, NPixY
%     p (used in ODE, but NOT here)
%     PixInptCtgrUse
%     C_SS_mean, C_CS_mean, C_IS_mean,
%     C_SC_mean, C_CC_mean, C_IC_mean,
%     C_SI_mean, C_CI_mean, C_II_mean

    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;

    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

    LDEIni.S = f(idxS);
    LDEIni.C = f(idxC);
    LDEIni.I = f(idxI);


    LDEUse = LDEIni;


    InhKillFlag = false; 

    % Default EKp / IKp (only used if you later set InhKillFlag = true)
    EKp.Thrsld   = 50;  EKp.Highist = [200,150]; EKp.HardBound = 200; EKp.Slope = 0;
    IKp.Thrsld   = 70;  IKp.Highist = [100,95];  IKp.HardBound = 120; IKp.Slope = 0.95;

    % Allow override from ParaODE if you ever want that:
    if isfield(ParaODE,'EKp'), EKp = ParaODE.EKp; end
    if isfield(ParaODE,'IKp'), IKp = ParaODE.IKp; end
    if isfield(ParaODE,'InhKillFlag'), InhKillFlag = ParaODE.InhKillFlag; end

    if InhKillFlag
        % From LDEIteration_16FuncMain_CombDom_0324
        LDEUse.E = LDEUse.S * (1-0.3077) + LDEUse.C * 0.3077;

        LDEUseEAdj = InhKill(LDEUse.E, ...
                             EKp.Thrsld, EKp.Highist, EKp.HardBound, EKp.Slope) ...
                     ./ LDEUse.E;
        LDEUse.S = LDEUse.S .* LDEUseEAdj;
        LDEUse.C = LDEUse.C .* LDEUseEAdj;

        LDEUse.I = InhKill(LDEUse.I, ...
                           IKp.Thrsld, IKp.Highist, IKp.HardBound, IKp.Slope);
        LDEUse.I = InhKill(LDEUse.I, ...
                           IKp.Thrsld, IKp.Highist, IKp.HardBound, IKp.Slope);
    end

    % ------------------------------------------------
    % 3. Compute "L4EUse" and "L4IUse" (synaptic input)
    %    EXACTLY as in LDEIteration_16FuncMain_CombDom_0324
    % ------------------------------------------------
    C_SS_mean = ParaODE.C_SS_mean;
    C_CS_mean = ParaODE.C_CS_mean;
    C_IS_mean = ParaODE.C_IS_mean;
    C_SC_mean = ParaODE.C_SC_mean;
    C_CC_mean = ParaODE.C_CC_mean;
    C_IC_mean = ParaODE.C_IC_mean;
    C_SI_mean = ParaODE.C_SI_mean;
    C_CI_mean = ParaODE.C_CI_mean;
    C_II_mean = ParaODE.C_II_mean;

    % Note: these lines are copied as-is from your LDEIteration code
    L4EUse = C_SS_mean*LDEUse.S + ...
             C_CS_mean*LDEUse.S + ...
             C_IS_mean*LDEUse.S + ...
             C_SC_mean*LDEUse.C + ...
             C_CC_mean*LDEUse.C + ...
             C_IC_mean*LDEUse.C;

    L4IUse = C_SI_mean*LDEUse.I + ...
             C_CI_mean*LDEUse.I + ...
             C_II_mean*LDEUse.I;

    % ------------------------------------------------
    % 4. Call the 16-function library once, with nan-check
    %    (this is your "instantaneous" 3D response)
    % ------------------------------------------------
    LDEOut = struct('S',[],'C',[],'I',[]);

    nanFlag = true;
    FuncN   = 3;   % as in original
    FuncUse = 1;   % start with finest

    while nanFlag && FuncUse <= FuncN
        LDEOut.S = LDEIterFunc_Grating_16Func_noFor_0324('S', ...
                         L4EUse, L4IUse, ...
                         ParaODE.PixInptCtgrUse);

        LDEOut.C = LDEIterFunc_Grating_16Func_noFor_0324('C', ...
                         L4EUse, L4IUse, ...
                         ParaODE.PixInptCtgrUse);

        LDEOut.I = LDEIterFunc_Grating_16Func_noFor_0324('I', ...
                         L4EUse, L4IUse, ...
                         ParaODE.PixInptCtgrUse);

        % Check for NaNs exactly like you intended
        nanFlag = any(isnan(LDEOut.S(:))) || ...
                  any(isnan(LDEOut.C(:))) || ...
                  any(isnan(LDEOut.I(:)));

        FuncUse = FuncUse + 1;
    end

    if nanFlag
        warning('Phi_from_2D_sigmoid: All function options produced NaN; returning input f as Φ(f).');
        Phi_f = f;
        return;
    end


    Phi_f = [LDEOut.S; LDEOut.C; LDEOut.I];
end







% function Phi_f = Phi_from_2D_sigmoid(f, ParaODE)
%     PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
%     idxS = 1:PixNumOut;
%     idxC = PixNumOut + (1:PixNumOut);
%     idxI = 2*PixNumOut + (1:PixNumOut);
% 
%     LDEIni.S = f(idxS);
%     LDEIni.C = f(idxC);
%     LDEIni.I = f(idxI);
% 
%     % Run exactly ONE epoch of the iteration from the current state f.
%     EpocOne = 1;
%     [f_output,~,~,~,~,~] = LDEIteration_16FuncMain_CombDom_0324( ...
%         ParaODE.PixInptCtgrUse, LDEIni, ParaODE.p, EpocOne, ...
%         ParaODE.C_SS_mean, ParaODE.C_CS_mean, ParaODE.C_IS_mean, ...
%         ParaODE.C_SC_mean, ParaODE.C_CC_mean, ParaODE.C_IC_mean, ...
%         ParaODE.C_SI_mean, ParaODE.C_CI_mean, ParaODE.C_II_mean, ...
%         ParaODE.N_HCOut, ParaODE.NPixX, ParaODE.NPixY, false, 'xn');
% 
%     % Traj1{1} = f_n, Traj1{2} = f_{n+1} (after one epoch with mixing p)
%     f_next = [f_output{2}.S; f_output{2}.C; f_output{2}.I];
% 
%     Phi_f = (f_next - (1 - ParaODE.p) * f) / ParaODE.p;
% end
