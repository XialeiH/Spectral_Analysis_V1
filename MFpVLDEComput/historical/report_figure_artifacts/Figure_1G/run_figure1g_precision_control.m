function run_figure1g_precision_control(sourceRoot,outputRoot,fieldHC,repeatId)
% Compare production and precision-matched CG/DNN response evaluation.
fieldHC=str2double(string(fieldHC)); repeatId=str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin'); addpath(fullfile(sourceRoot,'code'),'-begin'); addpath(fullfile(sourceRoot,'fast'),'-begin');
threads=str2double(getenv('SLURM_CPUS_PER_TASK')); if isfinite(threads)&&threads>=1,maxNumCompThreads(threads);end
B=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
R=load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat'));
fastSingle=h96_nn_set_condition_staticfirst(h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat')),100,0);
fastDouble=h96_nn_set_condition_staticfirst(h96_nn_prepare_base(B.ctx,'double',fullfile(sourceRoot,'fast','h96_models_double.mat')),100,0);
cgDouble=figure1g_prepare_cg_interpolants(R); cgSingle=prepare_single(R);
[recurrentE,recurrentI,l6]=initial_queries(B.initialState,B.ctx,R);
methods=["CG-Double","CG-Single","DNN-Double","DNN-Single"];
fns={@() cg_eval(recurrentE,recurrentI,l6,B.ctx.PixLGNCtgr,cgDouble,false), ...
     @() cg_eval(recurrentE,recurrentI,l6,B.ctx.PixLGNCtgr,cgSingle,true), ...
     @() h96_nn_all_population_responses_staticfirst(recurrentE,recurrentI,l6,fastDouble), ...
     @() h96_nn_all_population_responses_staticfirst(recurrentE,recurrentI,l6,fastSingle)};
reference=fns{1}(); rows=cell(4,1);
for methodId=1:4
    value=fns{methodId}(); [med,avg,sd]=timing(fns{methodId});
    rows{methodId}=table(fieldHC,repeatId,methods(methodId),med,avg,sd, ...
        norm(double(value)-reference,'fro')/max(norm(reference,'fro'),eps), ...
        max(abs(double(value(:))-reference(:))), ...
        'VariableNames',{'FieldHC','Repeat','Method','MedianSeconds','MeanSeconds','StdSeconds','RelativeError','MaxAbsoluteError'});
end
T=vertcat(rows{:}); out=fullfile(outputRoot,'precision',sprintf('precision_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function [recurrentE,recurrentI,l6]=initial_queries(state,ctx,R)
S=state.S;C=state.C;I=state.I;E=.6923*S+.3077*C;scale=L6Convert(E,ctx.EKpUse)./E;S=S.*scale;C=C.*scale;I=InhMulp(I,ctx.IKpUse);
recurrentE=[(ctx.C_SS_meanU*S+ctx.C_SC_meanU*C)/ctx.L4SEp;(ctx.C_CS_meanU*S+ctx.C_CC_meanU*C)/ctx.L4CEp;(ctx.C_IS_mean*S+ctx.C_IC_mean*C)/ctx.L4IEp];
recurrentI=[(ctx.C_SI_mean*I)/ctx.L4SIp;(ctx.C_CI_mean*I)/ctx.L4CIp;(ctx.C_II_mean*I)/ctx.L4IIp];
l6=repmat(median(R.l6Mesh),numel(S),1);
end

function cg=prepare_single(R)
l4E=single(unique(R.l4EMesh));l4I=single(unique(R.l4IMesh));l6=single(R.l6Mesh(:));names={'S','C','I'};cg.F=cell(5,3);
for p=1:3
    values=R.func.(names{p});
    for category=1:5,cg.F{category,p}=griddedInterpolant({l4I,l4E,l6},single(squeeze(values(category,:,:,:))),'linear','none');end
end
end

function output=cg_eval(recurrentE,recurrentI,l6,weights,cg,useSingle)
n=numel(l6);e=reshape(recurrentE,n,3);i=reshape(recurrentI,n,3);output=zeros(n,3);
if useSingle,e=single(e);i=single(i);l6=single(l6);weights=single(weights);end
for p=1:3,output(:,p)=figure1g_cg_population_response(e(:,p),i(:,p),l6,weights,cg.F(:,p));end
end

function [med,avg,sd]=timing(fn)
fn();samples=zeros(5,1);for k=1:5,timer=tic;fn();samples(k)=toc(timer);end;med=median(samples);avg=mean(samples);sd=std(samples);
end
