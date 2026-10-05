function run_figure1g_component_profile_full(sourceRoot,outputRoot,fieldHC,repeatId)
% Profile the shared and response-specific stages of CG and h96 DNN iterations.
fieldHC = str2double(string(fieldHC));
repeatId = str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin');
addpath(fullfile(sourceRoot,'code'),'-begin');
addpath(fullfile(sourceRoot,'fast'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1, maxNumCompThreads(threads); end

B = load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
    'ctx','initialState');
R = load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat'));
base = h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
fast = h96_nn_set_condition_symmetry(base,100,0);

cg = profile_cg(B.initialState,B.ctx,R,5);
dnn = profile_dnn(B.initialState,fast,5);
T = [profile_row(fieldHC,repeatId,"CG-Paper3",cg); ...
     profile_row(fieldHC,repeatId,"DNN-Class",dnn)];
out = fullfile(outputRoot,'component_profile', ...
    sprintf('profile_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function T = profile_row(fieldHC,repeatId,method,t)
responseFraction = t.Response/t.Total;
if method=="CG-Paper3"
    maxPossibleSpeedup = 1/max(1-responseFraction,eps);
else
    maxPossibleSpeedup = NaN;
end
T = table(fieldHC,repeatId,method,t.Total,t.Saturation,t.Recurrent,t.L6, ...
    t.Response,t.Scatter,t.Update,t.Allocation,responseFraction,maxPossibleSpeedup, ...
    'VariableNames',{'FieldHC','Repeat','Method','TotalSeconds','SaturationSeconds', ...
    'RecurrentSeconds','L6Seconds','ResponseSeconds','ScatterSeconds', ...
    'UpdateSeconds','AllocationConversionSeconds','ResponseFraction','MaxPossibleSpeedup'});
end

function t = profile_cg(initialState,ctx,R,iterations)
S=initialState.S; C=initialState.C; I=initialState.I;
fieldRows=ctx.N_HCOutY*ctx.NPixY; pars=ctx.L6pars{ctx.L6parId};
t=blank_timing(); total=tic;
for epoch=1:iterations
    q=tic; E=.6923*S+.3077*C; adjusted=L6Convert(E,ctx.EKpUse); scale=adjusted./E;
    Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,ctx.IKpUse); Euse=.6923*Suse+.3077*Cuse;
    t.Saturation=t.Saturation+toc(q);
    q=tic;
    rES=(ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
    rEC=(ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
    rEI=(ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp;
    rIS=(ctx.C_SI_mean*Iuse)/ctx.L4SIp; rIC=(ctx.C_CI_mean*Iuse)/ctx.L4CIp;
    rII=(ctx.C_II_mean*Iuse)/ctx.L4IIp; t.Recurrent=t.Recurrent+toc(q);
    q=tic; field=reshape(Euse,fieldRows,[]); padded=padarray(field,[1 1],'circular');
    filtered=conv2(padded,ctx.L6Kernel,'same'); l6=L6Convert(filtered(2:end-1,2:end-1),pars);
    l6=min(max(l6(:),3),3*numel(R.l6Mesh))./3; t.L6=t.L6+toc(q);
    q=tic;
    Sout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.S,rES,rIS,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    Cout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.C,rEC,rIC,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    Iout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.I,rEI,rII,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    t.Response=t.Response+toc(q);
    q=tic; S=ctx.p*Sout+(1-ctx.p)*S; C=ctx.p*Cout+(1-ctx.p)*C; I=ctx.p*Iout+(1-ctx.p)*I;
    t.Update=t.Update+toc(q);
end
t.Total=toc(total); t.Allocation=max(0,t.Total-t.Saturation-t.Recurrent-t.L6-t.Response-t.Update);
end

function t = profile_dnn(initialState,fast,iterations)
S=initialState.S; C=initialState.C; I=initialState.I; t=blank_timing(); total=tic;
for epoch=1:iterations
    q=tic; E=.6923*S+.3077*C; adjusted=L6Convert(E,fast.EKp); scale=adjusted./E;
    Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp); Euse=.6923*Suse+.3077*Cuse;
    t.Saturation=t.Saturation+toc(q);
    q=tic; recurrentE=fast.AE*[Suse;Cuse]; recurrentI=fast.AI*Iuse;
    t.Recurrent=t.Recurrent+toc(q);
    q=tic; field=reshape(Euse,fast.fieldRows,fast.fieldCols); padded=padarray(field,[1 1],'circular');
    filtered=conv2(padded,fast.kernel,'same'); l6Input=filtered(2:end-1,2:end-1);
    if isempty(fast.l6PP), l6=L6Convert(l6Input,fast.l6Pars); else, l6=ppval(fast.l6PP,l6Input); end
    l6=min(max(l6(:),3),fast.l6Max)./3; t.L6=t.L6+toc(q);
    q=tic; output=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fast);
    t.Response=t.Response+toc(q);
    q=tic; S=fast.p*output(:,1)+fast.oneMinusP*S; C=fast.p*output(:,2)+fast.oneMinusP*C;
    I=fast.p*output(:,3)+fast.oneMinusP*I; t.Update=t.Update+toc(q);
end
t.Total=toc(total); t.Allocation=max(0,t.Total-t.Saturation-t.Recurrent-t.L6-t.Response-t.Update);
end

function t=blank_timing()
t=struct('Total',0,'Saturation',0,'Recurrent',0,'L6',0,'Response',0,'Scatter',0,'Update',0,'Allocation',0);
end
