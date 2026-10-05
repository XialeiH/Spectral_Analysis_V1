function validate_figure1g_optimized_kernels(runRoot)
% Validate optimized kernels, then smoke-test field-size scaling.
repoRoot = fullfile(runRoot,'repo');
addpath(fullfile(repoRoot,'Utils'),'-begin');
addpath(fullfile(runRoot,'fast'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1
    maxNumCompThreads(threads);
end

R = load(fullfile(runRoot,'bundles','cg_response_c100_a0.mat'));
cgFast = figure1g_prepare_cg_interpolants(R);
modelFile = fullfile(runRoot,'fast','h96_models_double.mat');
fields = [4 20 40];
fprintf('VALIDATING OPTIMIZED FIGURE 1G KERNELS\n');
for fieldHC = fields
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
        'ctx','initialState');

    baseDouble = h96_nn_prepare_base(B.ctx,'double',modelFile);
    fastOld = h96_nn_set_condition(baseDouble,100,0);
    fastDouble = h96_nn_set_condition_staticfirst(baseDouble,100,0);
    [stateOld,oldInfo] = h96_nn_fixedpoint_compact(B.initialState,fastOld,3);
    [stateDouble,doubleInfo] = h96_nn_fixedpoint_staticfirst(B.initialState,fastDouble,3);
    dnnDoubleError = max(abs(pack(stateOld)-pack(stateDouble)));

    baseSingle = h96_nn_prepare_base(B.ctx,'single',modelFile);
    fastSingle = h96_nn_set_condition_staticfirst(baseSingle,100,0);
    [stateSingle,singleInfo] = h96_nn_fixedpoint_staticfirst(B.initialState,fastSingle,3);
    dnnSingleRelative = norm(pack(stateOld)-pack(stateSingle))/norm(pack(stateOld));
    [stateUnique,uniqueInfo] = h96_nn_fixedpoint_unique(B.initialState,fastSingle,3);
    dnnUniqueRelative = norm(pack(stateSingle)-pack(stateUnique))/norm(pack(stateSingle));
    fastSymmetry = h96_nn_set_condition_symmetry(baseSingle,100,0);
    [stateSymmetry,symmetryInfo] = h96_nn_fixedpoint_symmetry(B.initialState,fastSymmetry,3);
    dnnSymmetryRelative = norm(pack(stateSingle)-pack(stateSymmetry))/norm(pack(stateSingle));
    [stateAdaptive,adaptiveInfo] = h96_nn_fixedpoint_adaptive(B.initialState,fastSingle,3);
    dnnAdaptiveRelative = norm(pack(stateUnique)-pack(stateAdaptive))/norm(pack(stateUnique));

    [cgOld,cgOldSeconds] = cg_steps(B.initialState,B.ctx,R,[],1,false);
    [cgNew,cgNewSeconds] = cg_steps(B.initialState,B.ctx,R,cgFast,1,true);
    cgError = max(abs(pack(cgOld)-pack(cgNew)));

    fprintf(['FIELD=%d DNN_OLD=%.6f DNN_STATIC_DOUBLE=%.6f ', ...
        'DNN_STATIC_SINGLE=%.6f DNN_UNIQUE=%.6f DNN_SYMMETRY=%.6f ', ...
        'DNN_ADAPTIVE=%.6f DNN_DOUBLE_MAXERR=%.3e DNN_SINGLE_RELERR=%.3e DNN_UNIQUE_RELERR=%.3e ', ...
        'DNN_SYMMETRY_RELERR=%.3e UNIQUE_LAST=%d/%d/%d SYM_CLASSES=%d ', ...
        'DNN_ADAPTIVE_RELERR=%.3e ', ...
        'CG_OLD=%.6f CG_PREBUILT=%.6f CG_MAXERR=%.3e\n'], ...
        fieldHC,oldInfo.Seconds,doubleInfo.Seconds,singleInfo.Seconds, ...
        uniqueInfo.Seconds,symmetryInfo.Seconds,adaptiveInfo.Seconds,dnnDoubleError,dnnSingleRelative, ...
        dnnUniqueRelative,dnnSymmetryRelative, ...
        uniqueInfo.UniqueCounts(end,1),uniqueInfo.UniqueCounts(end,2), ...
        uniqueInfo.UniqueCounts(end,3),symmetryInfo.SymmetryClassCount, ...
        dnnAdaptiveRelative,cgOldSeconds,cgNewSeconds,cgError);
end
end

function [state,seconds] = cg_steps(initialState,ctx,R,cgFast,iterations,useFast)
S = initialState.S;
C = initialState.C;
I = initialState.I;
fieldRows = ctx.N_HCOutY*ctx.NPixY;
l6Pars = ctx.L6pars{ctx.L6parId};
timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S+0.3077.*C;
    adjustedE = L6Convert(E,ctx.EKpUse);
    scale = adjustedE./E;
    Suse = S.*scale;
    Cuse = C.*scale;
    Iuse = InhMulp(I,ctx.IKpUse);
    Euse = 0.6923.*Suse+0.3077.*Cuse;
    recurrentES = (ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
    recurrentEC = (ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
    recurrentEI = (ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp;
    recurrentIS = (ctx.C_SI_mean*Iuse)/ctx.L4SIp;
    recurrentIC = (ctx.C_CI_mean*Iuse)/ctx.L4CIp;
    recurrentII = (ctx.C_II_mean*Iuse)/ctx.L4IIp;
    field = reshape(Euse,fieldRows,[]);
    padded = padarray(field,[1 1],'circular');
    filtered = conv2(padded,ctx.L6Kernel,'same');
    l6Axis = L6Convert(filtered(2:end-1,2:end-1),l6Pars);
    l6Axis = min(max(l6Axis(:),3),3*numel(R.l6Mesh))./3;
    if useFast
        Sout = figure1g_cg_population_response(recurrentES,recurrentIS,l6Axis, ...
            ctx.PixLGNCtgr,cgFast.F(:,1));
        Cout = figure1g_cg_population_response(recurrentEC,recurrentIC,l6Axis, ...
            ctx.PixLGNCtgr,cgFast.F(:,2));
        Iout = figure1g_cg_population_response(recurrentEI,recurrentII,l6Axis, ...
            ctx.PixLGNCtgr,cgFast.F(:,3));
    else
        Sout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.S,recurrentES,recurrentIS,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
        Cout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.C,recurrentEC,recurrentIC,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
        Iout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
            R.func.I,recurrentEI,recurrentII,ctx.PixLGNCtgr,l6Axis,[1;2;3;4]);
    end
    S = ctx.p.*Sout+(1-ctx.p).*S;
    C = ctx.p.*Cout+(1-ctx.p).*C;
    I = ctx.p.*Iout+(1-ctx.p).*I;
end
seconds = toc(timer);
state = struct('S',S,'C',C,'I',I);
end

function x = pack(state)
x = [state.S;state.C;state.I];
end
