function D = Compute_D_QR_from_iteration_library(f, ParaODE, L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll)
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
    N_HCOut = ParaODE.N_HCOut;
    NPixX = ParaODE.NPixX;
    NPixY =ParaODE.NPixY;


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

   
    h = 1e-3;
    p = 0.33;
    EpocTest = 1;
  
    D = zeros(3*N, 3*N);

% Use Parallel Loop (Requires Parallel Computing Toolbox)
    for j = 1:3*N
    
        f_plus = f;
        f_minus = f;
    
        f_plus(j) = f(j) + h;
        f_minus(j) = f(j) - h;
    
        S_plus = f_plus(1:N);
        C_plus = f_plus(N+1:2*N);
        I_plus = f_plus(2*N+1:3*N);
    
        Ini_plus = struct('S', S_plus, 'C', C_plus, 'I', I_plus);

        S_minus = f_minus(1:N);
        C_minus = f_minus(N+1:2*N);
        I_minus = f_minus(2*N+1:3*N);
    
        Ini_minus = struct('S', S_minus, 'C', C_minus, 'I', I_minus);


            % iteration using precomputed library
        [LDEOutLib_plus,~,~,~,~,~] = ...
            LDEIteration_16FuncMain_CombDom_Phi( ...
            PixInptCtgrUse, Ini_plus, p, EpocTest, ...
            C_SS_mean, C_CS_mean, C_IS_mean, ...
            C_SC_mean, C_CC_mean, C_IC_mean, ...
            C_SI_mean, C_CI_mean, C_II_mean, ...
            L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, ...
            N_HCOut, NPixX, NPixY, false, 'xn');

        [LDEOutLib_minus,~,~,~,~,~] = ...
            LDEIteration_16FuncMain_CombDom_Phi( ...
            PixInptCtgrUse, Ini_minus, p, EpocTest, ...
            C_SS_mean, C_CS_mean, C_IS_mean, ...
            C_SC_mean, C_CC_mean, C_IC_mean, ...
            C_SI_mean, C_CI_mean, C_II_mean, ...
            L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, ...
            N_HCOut, NPixX, NPixY, false, 'xn');    


        Plus = [LDEOutLib_plus{end}.S; LDEOutLib_plus{end}.C; LDEOutLib_plus{end}.I];
        Minus = [LDEOutLib_minus{end}.S; LDEOutLib_minus{end}.C; LDEOutLib_minus{end}.I];



        Diff = (Plus - Minus) / (2*h);
    
    
        D(:, j) = Diff;
    
    end



end


