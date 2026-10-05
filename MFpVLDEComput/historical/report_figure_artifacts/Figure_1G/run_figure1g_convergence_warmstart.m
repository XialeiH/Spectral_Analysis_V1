function run_figure1g_convergence_warmstart(sourceRoot,outputRoot,fieldHC,repeatId)
% Compare cold CG/DNN solves with h96 initialization plus exact CG correction.
fieldHC=str2double(string(fieldHC)); repeatId=str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin');
addpath(fullfile(sourceRoot,'code'),'-begin'); addpath(fullfile(sourceRoot,'fast'),'-begin');
threads=str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads)&&threads>=1, maxNumCompThreads(threads); end
B=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
R=load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat'));
base=h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
fast=h96_nn_set_condition_symmetry(base,100,0);
tol=1e-6; consecutive=3; maxIterations=200;
initialState=B.initialState;
if repeatId>3
    rng(1729+repeatId,'twister');
    initialState.S=max(0,initialState.S.*exp(0.08*randn(size(initialState.S))-0.5*0.08^2));
    initialState.C=max(0,initialState.C.*exp(0.08*randn(size(initialState.C))-0.5*0.08^2));
    initialState.I=max(0,initialState.I.*exp(0.08*randn(size(initialState.I))-0.5*0.08^2));
end

[cgState,cgIter,cgTime,cgConv]=solve_state(initialState,@(x) step_cg(x,B.ctx,R),tol,consecutive,maxIterations);
[dnnState,dnnIter,dnnTime,dnnConv]=solve_state(initialState,@(x) step_dnn(x,fast),tol,consecutive,maxIterations);
[correctedState,correctionIter,correctionTime,correctionConv]=solve_state(dnnState,@(x) step_cg(x,B.ctx,R),tol,consecutive,maxIterations);
warmTotal=dnnTime+correctionTime;

rows=[result_row(fieldHC,repeatId,"CG-Cold",cgState,cgState,B.ctx,R,cgIter,cgTime,cgConv,0); ...
      result_row(fieldHC,repeatId,"DNN-Cold",dnnState,cgState,B.ctx,R,dnnIter,dnnTime,dnnConv,0); ...
      result_row(fieldHC,repeatId,"DNN-Init-CG-Correction",correctedState,cgState,B.ctx,R, ...
      correctionIter,warmTotal,correctionConv,dnnTime)];
out=fullfile(outputRoot,'convergence',sprintf('convergence_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(rows,out,'FileType','text','Delimiter','\t');
end

function [state,iterations,seconds,converged]=solve_state(initialState,stepFn,tol,needed,maxIterations)
state=initialState; stableCount=0; converged=false; timer=tic;
for iterations=1:maxIterations
    next=stepFn(state); residual=state_distance(next,state); state=next;
    if residual<tol, stableCount=stableCount+1; else, stableCount=0; end
    if stableCount>=needed, converged=true; break; end
end
seconds=toc(timer);
end

function T=result_row(fieldHC,repeatId,method,state,reference,ctx,R,iterations,seconds,converged,initializerSeconds)
next=step_cg(state,ctx,R); originalResidual=state_distance(next,state);
relativeStateError=state_distance(state,reference);
maxStateError=max(abs([state.S-reference.S;state.C-reference.C;state.I-reference.I]));
populationSError=norm(state.S-reference.S)/max(norm(reference.S),eps);
populationCError=norm(state.C-reference.C)/max(norm(reference.C),eps);
populationIError=norm(state.I-reference.I)/max(norm(reference.I),eps);
T=table(fieldHC,repeatId,method,iterations,seconds,converged,initializerSeconds, ...
    originalResidual,relativeStateError,maxStateError,populationSError,populationCError,populationIError, ...
    max(state.S),max(state.C),max(state.I), ...
    'VariableNames',{'FieldHC','Repeat','Method','Iterations','Seconds','Converged', ...
    'InitializerSeconds','OriginalCGResidual','RelativeStateError','MaxStateError', ...
    'SError','CError','IError','MaxS','MaxC','MaxI'});
end

function d=state_distance(a,b)
delta=[a.S-b.S;a.C-b.C;a.I-b.I]; base=[b.S;b.C;b.I];
d=norm(delta)/max(norm(base),eps);
end

function next=step_cg(state,ctx,R)
S=state.S; C=state.C; I=state.I; fieldRows=ctx.N_HCOutY*ctx.NPixY;
E=.6923*S+.3077*C; scale=L6Convert(E,ctx.EKpUse)./E;
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,ctx.IKpUse); Euse=.6923*Suse+.3077*Cuse;
rES=(ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
rEC=(ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
rEI=(ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp;
rIS=(ctx.C_SI_mean*Iuse)/ctx.L4SIp; rIC=(ctx.C_CI_mean*Iuse)/ctx.L4CIp; rII=(ctx.C_II_mean*Iuse)/ctx.L4IIp;
field=reshape(Euse,fieldRows,[]); filtered=conv2(padarray(field,[1 1],'circular'),ctx.L6Kernel,'same');
l6=L6Convert(filtered(2:end-1,2:end-1),ctx.L6pars{ctx.L6parId}); l6=min(max(l6(:),3),3*numel(R.l6Mesh))./3;
Sout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.S,rES,rIS,ctx.PixLGNCtgr,l6,[1;2;3;4]);
Cout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.C,rEC,rIC,ctx.PixLGNCtgr,l6,[1;2;3;4]);
Iout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.I,rEI,rII,ctx.PixLGNCtgr,l6,[1;2;3;4]);
next=struct('S',ctx.p*Sout+(1-ctx.p)*S,'C',ctx.p*Cout+(1-ctx.p)*C,'I',ctx.p*Iout+(1-ctx.p)*I);
end

function next=step_dnn(state,fast)
S=state.S; C=state.C; I=state.I;
E=.6923*S+.3077*C; scale=L6Convert(E,fast.EKp)./E;
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp); Euse=.6923*Suse+.3077*Cuse;
recurrentE=fast.AE*[Suse;Cuse]; recurrentI=fast.AI*Iuse;
field=reshape(Euse,fast.fieldRows,fast.fieldCols); filtered=conv2(padarray(field,[1 1],'circular'),fast.kernel,'same');
l6Input=filtered(2:end-1,2:end-1);
if isempty(fast.l6PP), l6=L6Convert(l6Input,fast.l6Pars); else, l6=ppval(fast.l6PP,l6Input); end
l6=min(max(l6(:),3),fast.l6Max)./3;
output=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fast);
next=struct('S',fast.p*output(:,1)+fast.oneMinusP*S, ...
    'C',fast.p*output(:,2)+fast.oneMinusP*C,'I',fast.p*output(:,3)+fast.oneMinusP*I);
end
