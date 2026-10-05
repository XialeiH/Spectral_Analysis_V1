function profile_figure1g_cg_iteration(runRoot)
% Measure the original Paper 3 operations within a CG iteration.
addpath(fullfile(runRoot,'repo','Utils'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1
    maxNumCompThreads(threads);
end
R = load(fullfile(runRoot,'bundles','cg_response_c100_a0.mat'));
for fieldHC = [4 40]
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
        'ctx','initialState');
    ctx = B.ctx;
    S = B.initialState.S; C = B.initialState.C; I = B.initialState.I;
    fieldRows = ctx.N_HCOutY*ctx.NPixY;
    l6Pars = ctx.L6pars{ctx.L6parId};
    iterations = 5;
    timing = zeros(1,5);
    totalTimer = tic;
    for epoch = 1:iterations
        timer = tic;
        E = 0.6923.*S+0.3077.*C;
        adjustedE = L6Convert(E,ctx.EKpUse);
        scale = adjustedE./E;
        Suse = S.*scale; Cuse = C.*scale; Iuse = InhMulp(I,ctx.IKpUse);
        Euse = 0.6923.*Suse+0.3077.*Cuse;
        timing(1) = timing(1)+toc(timer);

        timer = tic;
        recurrentES = (ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
        recurrentEC = (ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
        recurrentEI = (ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp;
        recurrentIS = (ctx.C_SI_mean*Iuse)/ctx.L4SIp;
        recurrentIC = (ctx.C_CI_mean*Iuse)/ctx.L4CIp;
        recurrentII = (ctx.C_II_mean*Iuse)/ctx.L4IIp;
        timing(2) = timing(2)+toc(timer);

        timer = tic;
        field = reshape(Euse,fieldRows,[]);
        padded = padarray(field,[1 1],'circular');
        filtered = conv2(padded,ctx.L6Kernel,'same');
        l6Axis = L6Convert(filtered(2:end-1,2:end-1),l6Pars);
        l6Axis = min(max(l6Axis(:),3),3*numel(R.l6Mesh))./3;
        timing(3) = timing(3)+toc(timer);

        timer = tic;
        Sout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.S,recurrentES,recurrentIS,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
        Cout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.C,recurrentEC,recurrentIC,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
        Iout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.I,recurrentEI,recurrentII,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
        timing(4) = timing(4)+toc(timer);

        timer = tic;
        S = ctx.p.*Sout+(1-ctx.p).*S;
        C = ctx.p.*Cout+(1-ctx.p).*C;
        I = ctx.p.*Iout+(1-ctx.p).*I;
        timing(5) = timing(5)+toc(timer);
    end
    total = toc(totalTimer);
    fprintf(['FIELD=%d PER_ITER_TOTAL=%.6f SATURATION=%.6f MATRIX=%.6f ', ...
        'L6=%.6f LOOKUP=%.6f UPDATE=%.6f LOOKUP_PCT=%.2f MATRIX_PCT=%.2f\n'], ...
        fieldHC,total/iterations,timing(1)/iterations,timing(2)/iterations, ...
        timing(3)/iterations,timing(4)/iterations,timing(5)/iterations, ...
        100*timing(4)/total,100*timing(2)/total);
end
end
