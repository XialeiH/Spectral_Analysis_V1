function run_figure1g_field_trial(runRoot, fieldHC, repeatId)
% Time paired 50-iteration Paper 3 CG and optimized h96 DNN solves.
fieldHC = str2double(string(fieldHC));
repeatId = str2double(string(repeatId));
iterations = 50;
repoRoot = fullfile(runRoot, 'repo');
addpath(fullfile(repoRoot,'Utils'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
addpath(fullfile(runRoot,'fast'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads >= 1
    maxNumCompThreads(threads);
end

B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
    'ctx','initialState');
persistent cachedRunRoot responseCache fastByField
if isempty(cachedRunRoot) || ~strcmp(cachedRunRoot,runRoot)
    responseCache = load(fullfile(runRoot,'bundles','cg_response_c100_a0.mat'));
    fastByField = containers.Map('KeyType','double','ValueType','any');
    cachedRunRoot = runRoot;
end
R = responseCache;
if isKey(fastByField,fieldHC)
    fast = fastByField(fieldHC);
else
    modelFile = fullfile(runRoot,'fast','h96_models_double.mat');
    base = h96_nn_prepare_base(B.ctx,'single',modelFile);
    fast = h96_nn_set_condition_symmetry(base,100,0);
    fastByField(fieldHC) = fast;
end

if mod(repeatId,2)==1
    cg = run_cg(B.initialState,B.ctx,R,iterations);
    dnn = run_dnn(B.initialState,fast,iterations);
else
    dnn = run_dnn(B.initialState,fast,iterations);
    cg = run_cg(B.initialState,B.ctx,R,iterations);
end

T = table(fieldHC,repeatId,iterations,cg.Seconds,dnn.Seconds, ...
    cg.Seconds/iterations,dnn.Seconds/iterations,cg.Checksum,dnn.Checksum, ...
    'VariableNames',{'FieldHC','Repeat','Iterations','CGSeconds','DNNSeconds', ...
    'CGSecondsPerIteration','DNNSecondsPerIteration','CGChecksum','DNNChecksum'});
out = fullfile(runRoot,'results',sprintf('field_%02d_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
fprintf('FIELD=%d REPEAT=%d CG=%.6f DNN=%.6f SPEEDUP=%.3f\n', ...
    fieldHC,repeatId,cg.Seconds,dnn.Seconds,cg.Seconds/dnn.Seconds);
end

function info = run_dnn(initialState,fast,iterations)
timer = tic;
[state,~] = h96_nn_fixedpoint_symmetry(initialState,fast,iterations);
info.Seconds = toc(timer);
info.Checksum = sum(state.S)+sum(state.C)+sum(state.I);
end

function info = run_cg(initialState,ctx,R,iterations)
S = initialState.S; C = initialState.C; I = initialState.I;
fieldRows = ctx.N_HCOutY*ctx.NPixY;
l6Pars = ctx.L6pars{ctx.L6parId};
categoryOrder = [1;2;3;4];
timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S + 0.3077.*C;
    adjustedE = L6Convert(E,ctx.EKpUse);
    adjustment = adjustedE./E;
    Suse = S.*adjustment;
    Cuse = C.*adjustment;
    Iuse = InhMulp(I,ctx.IKpUse);
    Euse = 0.6923.*Suse + 0.3077.*Cuse;
    recurrentES = (ctx.C_SS_meanU*Suse + ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
    recurrentEC = (ctx.C_CS_meanU*Suse + ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
    recurrentEI = (ctx.C_IS_mean*Suse + ctx.C_IC_mean*Cuse)/ctx.L4IEp;
    recurrentIS = (ctx.C_SI_mean*Iuse)/ctx.L4SIp;
    recurrentIC = (ctx.C_CI_mean*Iuse)/ctx.L4CIp;
    recurrentII = (ctx.C_II_mean*Iuse)/ctx.L4IIp;
    field = reshape(Euse,fieldRows,[]);
    padded = padarray(field,[1 1],'circular');
    filtered = conv2(padded,ctx.L6Kernel,'same');
    l6Index = L6Convert(filtered(2:end-1,2:end-1),l6Pars);
    l6Index = min(max(l6Index(:),3),3*numel(R.l6Mesh))./3;
    Sout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
        R.func.S,recurrentES,recurrentIS,ctx.PixLGNCtgr,l6Index,categoryOrder);
    Cout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
        R.func.C,recurrentEC,recurrentIC,ctx.PixLGNCtgr,l6Index,categoryOrder);
    Iout = LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
        R.func.I,recurrentEI,recurrentII,ctx.PixLGNCtgr,l6Index,categoryOrder);
    S = ctx.p.*Sout + (1-ctx.p).*S;
    C = ctx.p.*Cout + (1-ctx.p).*C;
    I = ctx.p.*Iout + (1-ctx.p).*I;
end
info.Seconds = toc(timer);
info.Checksum = sum(S)+sum(C)+sum(I);
end
