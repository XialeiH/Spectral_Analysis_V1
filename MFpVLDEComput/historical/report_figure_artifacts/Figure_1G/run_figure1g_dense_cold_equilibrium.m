function run_figure1g_dense_cold_equilibrium(sourceRoot,outputRoot,fieldHC,repeatId)
% Paired cold-start equilibrium timing with dense recurrent multiplication.
fieldHC = str2double(string(fieldHC));
repeatId = str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin');
addpath(fullfile(sourceRoot,'code'),'-begin');
addpath(fullfile(sourceRoot,'fast'),'-begin');
addpath(fullfile(outputRoot,'code'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1, maxNumCompThreads(threads); end

B = load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
    'ctx','initialState');
R = load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat'));
base = h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
fast = h96_nn_set_condition_symmetry(base,100,0);
[~,pairColumn,pairWeight] = find(fast.mixAggregate);
fast.mixWeight = zeros(fast.activePairCount,1);
fast.mixWeight(pairColumn) = pairWeight;
fast = rmfield(fast,'mixAggregate');

tol = 1e-6;
consecutive = 3;
maxIterations = 200;
initialState = B.initialState;
if repeatId>3
    rng(1729+repeatId,'twister');
    initialState.S = max(0,initialState.S.*exp(0.08*randn(size(initialState.S))-0.5*0.08^2));
    initialState.C = max(0,initialState.C.*exp(0.08*randn(size(initialState.C))-0.5*0.08^2));
    initialState.I = max(0,initialState.I.*exp(0.08*randn(size(initialState.I))-0.5*0.08^2));
end

if mod(repeatId,2)==1
    cg = solve_state(initialState,@(x) step_cg_dense(x,B.ctx,R),tol,consecutive,maxIterations);
    dnn = solve_state(initialState,@(x) step_dnn_dense(x,fast),tol,consecutive,maxIterations);
else
    dnn = solve_state(initialState,@(x) step_dnn_dense(x,fast),tol,consecutive,maxIterations);
    cg = solve_state(initialState,@(x) step_cg_dense(x,B.ctx,R),tol,consecutive,maxIterations);
end

T = [result_row(fieldHC,repeatId,"CG",cg); ...
     result_row(fieldHC,repeatId,"DNN surrogate",dnn)];
out = fullfile(outputRoot,'results', ...
    sprintf('dense_cold_%02dHC_repeat_%03d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
fprintf('DENSE_COLD field=%d repeat=%d CG=%.6f(%d) DNN=%.6f(%d)\n', ...
    fieldHC,repeatId,cg.Seconds,cg.Iterations,dnn.Seconds,dnn.Iterations);
end

function info = solve_state(initialState,stepFn,tol,needed,maxIterations)
state = initialState;
stableCount = 0;
converged = false;
timer = tic;
for iterations = 1:maxIterations
    next = stepFn(state);
    residual = state_distance(next,state);
    state = next;
    if residual<tol, stableCount=stableCount+1; else, stableCount=0; end
    if stableCount>=needed, converged=true; break; end
end
info.Seconds = toc(timer);
info.Iterations = iterations;
info.Converged = converged;
info.FinalResidual = residual;
info.Checksum = sum(state.S)+sum(state.C)+sum(state.I);
end

function T = result_row(fieldHC,repeatId,method,info)
T = table(fieldHC,repeatId,method,info.Iterations,info.Seconds, ...
    info.Converged,info.FinalResidual,info.Checksum, ...
    'VariableNames',{'FieldHC','Repeat','Method','Iterations','Seconds', ...
    'Converged','FinalResidual','Checksum'});
end

function next = step_cg_dense(state,ctx,R)
S=state.S; C=state.C; I=state.I;
fieldRows=ctx.N_HCOutY*ctx.NPixY;
E=.6923*S+.3077*C;
scale=L6Convert(E,ctx.EKpUse)./E;
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,ctx.IKpUse);
Euse=.6923*Suse+.3077*Cuse;
rES=(dense_matvec(ctx.C_SS_meanU,Suse)+dense_matvec(ctx.C_SC_meanU,Cuse))/ctx.L4SEp;
rEC=(dense_matvec(ctx.C_CS_meanU,Suse)+dense_matvec(ctx.C_CC_meanU,Cuse))/ctx.L4CEp;
rEI=(dense_matvec(ctx.C_IS_mean,Suse)+dense_matvec(ctx.C_IC_mean,Cuse))/ctx.L4IEp;
rIS=dense_matvec(ctx.C_SI_mean,Iuse)/ctx.L4SIp;
rIC=dense_matvec(ctx.C_CI_mean,Iuse)/ctx.L4CIp;
rII=dense_matvec(ctx.C_II_mean,Iuse)/ctx.L4IIp;
field=reshape(Euse,fieldRows,[]);
filtered=conv2(padarray(field,[1 1],'circular'),ctx.L6Kernel,'same');
l6=L6Convert(filtered(2:end-1,2:end-1),ctx.L6pars{ctx.L6parId});
l6=min(max(l6(:),3),3*numel(R.l6Mesh))./3;
Sout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
    R.func.S,rES,rIS,ctx.PixLGNCtgr,l6,[1;2;3;4]);
Cout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
    R.func.C,rEC,rIC,ctx.PixLGNCtgr,l6,[1;2;3;4]);
Iout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh, ...
    R.func.I,rEI,rII,ctx.PixLGNCtgr,l6,[1;2;3;4]);
next=struct('S',ctx.p*Sout+(1-ctx.p)*S, ...
    'C',ctx.p*Cout+(1-ctx.p)*C,'I',ctx.p*Iout+(1-ctx.p)*I);
end

function next = step_dnn_dense(state,fast)
S=state.S; C=state.C; I=state.I;
E=.6923*S+.3077*C;
scale=L6Convert(E,fast.EKp)./E;
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp);
Euse=.6923*Suse+.3077*Cuse;
recurrentE=dense_matvec(fast.AE,[Suse;Cuse]);
recurrentI=dense_matvec(fast.AI,Iuse);
field=reshape(Euse,fast.fieldRows,fast.fieldCols);
filtered=conv2(padarray(field,[1 1],'circular'),fast.kernel,'same');
l6Input=filtered(2:end-1,2:end-1);
if isempty(fast.l6PP), l6=L6Convert(l6Input,fast.l6Pars); else, l6=ppval(fast.l6PP,l6Input); end
l6=min(max(l6(:),3),fast.l6Max)./3;
output=h96_nn_all_population_responses_symmetry_nosparse(recurrentE,recurrentI,l6,fast);
next=struct('S',fast.p*output(:,1)+fast.oneMinusP*S, ...
    'C',fast.p*output(:,2)+fast.oneMinusP*C, ...
    'I',fast.p*output(:,3)+fast.oneMinusP*I);
end

function y = dense_matvec(A,x)
% Sparse storage is only for memory; every block multiplication is dense.
nRows=size(A,1); nCols=size(A,2);
targetBytes=768*1024^2;
rowsPerBlock=max(1,floor(targetBytes/(8*nCols)));
y=zeros(nRows,1);
for first=1:rowsPerBlock:nRows
    rows=first:min(first+rowsPerBlock-1,nRows);
    denseBlock=full(A(rows,:));
    y(rows)=denseBlock*x;
end
end

function d = state_distance(a,b)
delta=[a.S-b.S;a.C-b.C;a.I-b.I];
base=[b.S;b.C;b.I];
d=norm(delta)/max(norm(base),eps);
end
