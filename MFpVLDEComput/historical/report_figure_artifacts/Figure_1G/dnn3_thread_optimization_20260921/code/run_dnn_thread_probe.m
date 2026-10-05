function run_dnn_thread_probe(sourceRoot,outputRoot,repeatId,nnThreads,phase)
% Same frozen 3x3 DNN update; vary only MATLAB computational thread limit.
repeatId = str2double(string(repeatId));
nnThreads = str2double(string(nnThreads));
assert(ismember(nnThreads,[1 2 4 8 16]));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin');
addpath(fullfile(sourceRoot,'code'),'-begin');
addpath(fullfile(sourceRoot,'fast'),'-begin');
addpath(fullfile(outputRoot,'code'),'-begin');
maxNumCompThreads(nnThreads);
assert(maxNumCompThreads==nnThreads);
B = load(fullfile(sourceRoot,'bundles','field_03HC.mat'),'ctx','initialState');
initialState = B.initialState;
if repeatId>3
    rng(1729+repeatId,'twister');
    initialState.S = max(0,initialState.S.*exp(0.08*randn(size(initialState.S))-0.5*0.08^2));
    initialState.C = max(0,initialState.C.*exp(0.08*randn(size(initialState.C))-0.5*0.08^2));
    initialState.I = max(0,initialState.I.*exp(0.08*randn(size(initialState.I))-0.5*0.08^2));
end
base = h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
fast = h96_nn_set_condition_symmetry(base,100,0);
stepFn = @(x) step_dnn_fast(x,fast);
tol = 5e-3; needed = 3; maxIterations = 200;
state = initialState; stableCount = 0; converged = false;
timer = tic;
for iterations = 1:maxIterations
    next = stepFn(state);
    residual = state_distance(next,state);
    state = next;
    if residual<tol, stableCount=stableCount+1; else, stableCount=0; end
    if stableCount>=needed, converged=true; break; end
end
seconds = toc(timer);
checksum = sum(state.S)+sum(state.C)+sum(state.I);
allocatedCPUs = str2double(getenv('SLURM_CPUS_PER_TASK'));
node = string(getenv('SLURMD_NODENAME'));
T = table(3,repeatId,"DNN surrogate",iterations,seconds,converged,residual, ...
    checksum,nnThreads,allocatedCPUs,node, ...
    'VariableNames',{'FieldHC','Repeat','Method','Iterations','Seconds', ...
    'Converged','FinalResidual','Checksum','NNThreads','AllocatedCPUs','Node'});
stem = sprintf('dnn_03HC_threads_%02d_repeat_%03d',nnThreads,repeatId);
% Validation and output are outside the original iteration-only timer.
save(fullfile(outputRoot,char(phase),[stem '.mat']),'state','T');
writetable(T,fullfile(outputRoot,char(phase),[stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
fprintf('threads=%d repeat=%d seconds=%.9f iterations=%d residual=%.9g\n', ...
    nnThreads,repeatId,seconds,iterations,residual);
end

function next = step_dnn_fast(state,fast)
S=state.S; C=state.C; I=state.I;
E=.6923*S+.3077*C;
scale=L6Convert(E,fast.EKp)./E;
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp);
Euse=.6923*Suse+.3077*Cuse;
recurrentE=fast.AE*[Suse;Cuse];
recurrentI=fast.AI*Iuse;
field=reshape(Euse,fast.fieldRows,fast.fieldCols);
filtered=conv2(padarray(field,[1 1],'circular'),fast.kernel,'same');
l6Input=filtered(2:end-1,2:end-1);
if isempty(fast.l6PP), l6=L6Convert(l6Input,fast.l6Pars); else, l6=ppval(fast.l6PP,l6Input); end
l6=min(max(l6(:),3),fast.l6Max)./3;
output=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fast);
next=struct('S',fast.p*output(:,1)+fast.oneMinusP*S, ...
    'C',fast.p*output(:,2)+fast.oneMinusP*C, ...
    'I',fast.p*output(:,3)+fast.oneMinusP*I);
end

function d = state_distance(a,b)
delta=[a.S-b.S;a.C-b.C;a.I-b.I];
base=[b.S;b.C;b.I];
d=norm(delta)/max(norm(base),eps);
end
