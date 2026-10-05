function out = Compare_forward_direct_vs_library(f, ParaODE, L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, p, EpocTest)
%COMPARE_FORWARD_DIRECT_VS_LIBRARY
% Compare one-shot sigmoid forward (direct) vs. one-step library iteration.
%
% Inputs:
%   f              - 3N vector [S; C; I]
%   ParaODE        - struct with connectivity and PixInptCtgrUse
%   L4EmeshXAll    - mesh grid for library (as used in iteration library)
%   L4ImeshYAll    - mesh grid for library (as used in iteration library)
%   LDEFrfuncAll   - precomputed library response functions
%   p              - library mixing parameter (e.g., 0.33)
%   EpocTest       - number of library iterations (use 1 to match Jacobian)
%
% Output:
%   out.direct     - struct with fields S,C,I from direct forward
%   out.library    - struct with fields S,C,I from library forward
%   out.diff_norm  - norm of differences for each population
%   out.diff_max   - max abs difference for each population

    if nargin < 7
        EpocTest = 1;
    end

    % Basic sizes and indexing
    PixNumOut = ParaODE.N_HCOut^2 * ParaODE.NPixX * ParaODE.NPixY;
    N         = PixNumOut;
    idxS = 1:PixNumOut;
    idxC = PixNumOut + (1:PixNumOut);
    idxI = 2*PixNumOut + (1:PixNumOut);

    % Connectivity
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

    % Split f
    S = f(idxS);
    C = f(idxC);
    I = f(idxI);

    %% Direct forward (one-shot sigmoid)
    L4E_dir = C_SS_mean*S + C_CS_mean*S + C_IS_mean*S + ...
              C_SC_mean*C + C_CC_mean*C + C_IC_mean*C;
    L4I_dir = C_SI_mean*I + C_CI_mean*I + C_II_mean*I;

    Psi_S_dir = local_Psi('S', L4E_dir, L4I_dir, PixInptCtgrUse);
    Psi_C_dir = local_Psi('C', L4E_dir, L4I_dir, PixInptCtgrUse);
    Psi_I_dir = local_Psi('I', L4E_dir, L4I_dir, PixInptCtgrUse);

    %% Library forward (one iteration)
    Ini = struct('S', S, 'C', C, 'I', I);
    N_HCOut = ParaODE.N_HCOut;
    NPixX = ParaODE.NPixX;
    NPixY = ParaODE.NPixY;

    [LDEOutLib, ~, ~, ~, ~, ~] = LDEIteration_16FuncMain_CombDom_Phi( ...
        PixInptCtgrUse, Ini, p, EpocTest, ...
        C_SS_mean, C_CS_mean, C_IS_mean, ...
        C_SC_mean, C_CC_mean, C_IC_mean, ...
        C_SI_mean, C_CI_mean, C_II_mean, ...
        L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, ...
        N_HCOut, NPixX, NPixY, false, 'xn');
 
    Psi_S_lib = LDEOutLib{end}.S(:);
    Psi_C_lib = LDEOutLib{end}.C(:);
    Psi_I_lib = LDEOutLib{end}.I(:);

    %% Differences
    diff_S = Psi_S_dir - Psi_S_lib;
    diff_C = Psi_C_dir - Psi_C_lib;
    diff_I = Psi_I_dir - Psi_I_lib;

    out.direct.S = Psi_S_dir;
    out.direct.C = Psi_C_dir;
    out.direct.I = Psi_I_dir;

    out.library.S = Psi_S_lib;
    out.library.C = Psi_C_lib;
    out.library.I = Psi_I_lib;

    out.diff_norm.S = norm(diff_S);
    out.diff_norm.C = norm(diff_C);
    out.diff_norm.I = norm(diff_I);

    out.diff_max.S = max(abs(diff_S));
    out.diff_max.C = max(abs(diff_C));
    out.diff_max.I = max(abs(diff_I));

    fprintf('Direct vs Library (Epoc=%d, p=%.3f) diff norms: S=%g, C=%g, I=%g\n', ...
        EpocTest, p, out.diff_norm.S, out.diff_norm.C, out.diff_norm.I);
    fprintf('Direct vs Library (Epoc=%d, p=%.3f) max abs:   S=%g, C=%g, I=%g\n', ...
        EpocTest, p, out.diff_max.S, out.diff_max.C, out.diff_max.I);
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
