function inspect_field_scaling(runRoot,outputRoot)
% Side diagnostic only: retain dense CG matvecs and original interpolation.
% Balanced field order; profiling results do not replace benchmark datapoints.
addpath(fullfile(runRoot,'repo','Utils'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
addpath(fullfile(runRoot,'fast'),'-begin');
addpath(fullfile(outputRoot,'code'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1
    maxNumCompThreads(threads);
end
R = load(fullfile(runRoot,'bundles','cg_response_c100_a0.mat'));
fprintf('LIBRARY_GRID %d x %d x %d; CATEGORIES %d; LOOKUPS_PER_ITER %d\n', ...
    numel(unique(R.l4EMesh)),numel(unique(R.l4IMesh)),numel(R.l6Mesh), ...
    size(R.func.S,1),size(R.func.S,1)+size(R.func.C,1)+size(R.func.I,1));
for h = [3 4]
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',h)));
    fast = h96_nn_set_condition_symmetry(h96_nn_prepare_base(B.ctx,'single', ...
        fullfile(runRoot,'fast','h96_models_double.mat')),100,0);
    fprintf('STRUCTURE FIELD=%d PIXELS=%d ACTIVE_PAIRS=%d NN_CLASSES=%d AE_SPARSE=%d AI_SPARSE=%d\n', ...
        h,fast.n,fast.activePairCount,fast.symClassCount,issparse(fast.AE),issparse(fast.AI));
end
for fieldHC = [3 4 4 3]
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
        'ctx','initialState');
    ctx = B.ctx;
    S = B.initialState.S; C = B.initialState.C; I = B.initialState.I;
    fieldRows = ctx.N_HCOutY*ctx.NPixY;
    l6Pars = ctx.L6pars{ctx.L6parId};
    iterations = 10;
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
        recurrentES = (dense_matvec(ctx.C_SS_meanU,Suse)+dense_matvec(ctx.C_SC_meanU,Cuse))/ctx.L4SEp;
        recurrentEC = (dense_matvec(ctx.C_CS_meanU,Suse)+dense_matvec(ctx.C_CC_meanU,Cuse))/ctx.L4CEp;
        recurrentEI = (dense_matvec(ctx.C_IS_mean,Suse)+dense_matvec(ctx.C_IC_mean,Cuse))/ctx.L4IEp;
        recurrentIS = (dense_matvec(ctx.C_SI_mean,Iuse))/ctx.L4SIp;
        recurrentIC = (dense_matvec(ctx.C_CI_mean,Iuse))/ctx.L4CIp;
        recurrentII = (dense_matvec(ctx.C_II_mean,Iuse))/ctx.L4IIp;
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

function y = dense_matvec(A,x)
nRows=size(A,1); nCols=size(A,2);
rowsPerBlock=max(1,floor((768*1024^2)/(8*nCols)));
y=zeros(nRows,1);
for first=1:rowsPerBlock:nRows
    rows=first:min(first+rowsPerBlock-1,nRows);
    y(rows)=full(A(rows,:))*x;
end
end

