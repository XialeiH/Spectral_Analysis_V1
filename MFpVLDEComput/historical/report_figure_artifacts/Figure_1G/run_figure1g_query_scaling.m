function run_figure1g_query_scaling(sourceRoot,outputRoot,fieldHC,repeatId)
% Measure response cost as the number of distinct local queries increases.
fieldHC=str2double(string(fieldHC)); repeatId=str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin'); addpath(fullfile(sourceRoot,'code'),'-begin');
addpath(fullfile(sourceRoot,'fast'),'-begin');
threads=str2double(getenv('SLURM_CPUS_PER_TASK')); if isfinite(threads)&&threads>=1, maxNumCompThreads(threads); end
B=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
R=load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat'));
base=h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
fast=h96_nn_set_condition_symmetry(base,100,0); cg=figure1g_prepare_cg_interpolants(R);
q=record_query(B.initialState,B.ctx,R,25); n=numel(q.l6);
regimes=["Fixed","SqrtN","Full"];
requestedK=[min(400,n),max(1,ceil(sqrt(n))),n];
methods=["CG-Optimized-Pixel","CG-Optimized-Class","DNN-Pixel","DNN-Class"];
rows=cell(numel(regimes)*numel(methods),1); rowId=0;
for regimeId=1:numel(regimes)
    qr=make_queries(q,requestedK(regimeId),R);
    evaluators={@() cg_pixel(qr,B.ctx,cg),@() cg_class(qr,B.ctx,cg), ...
        @() h96_nn_all_population_responses_staticfirst(qr.recurrentE,qr.recurrentI,qr.l6,fast), ...
        @() h96_nn_all_population_responses_unique(qr.recurrentE,qr.recurrentI,qr.l6,fast)};
    for methodId=1:numel(methods)
        fn=evaluators{methodId}; fn(); samples=zeros(5,1);
        for k=1:numel(samples), timer=tic; fn(); samples(k)=toc(timer); end
        rowId=rowId+1;
        rows{rowId}=table(fieldHC,repeatId,regimes(regimeId),requestedK(regimeId), ...
            qr.actualK,methods(methodId),median(samples),mean(samples),std(samples),n, ...
            'VariableNames',{'FieldHC','Repeat','Regime','RequestedK','ActualK', ...
            'Method','MedianSeconds','MeanSeconds','StdSeconds','PixelCount'});
    end
end
T=vertcat(rows{:}); out=fullfile(outputRoot,'query_scaling',sprintf('query_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function q=record_query(initialState,ctx,R,targetEpoch)
S=initialState.S; C=initialState.C; I=initialState.I; fieldRows=ctx.N_HCOutY*ctx.NPixY;
for epoch=1:targetEpoch
    E=.6923*S+.3077*C; scale=L6Convert(E,ctx.EKpUse)./E; Suse=S.*scale; Cuse=C.*scale;
    Iuse=InhMulp(I,ctx.IKpUse); Euse=.6923*Suse+.3077*Cuse;
    rES=(ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp; rEC=(ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
    rEI=(ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp; rIS=(ctx.C_SI_mean*Iuse)/ctx.L4SIp;
    rIC=(ctx.C_CI_mean*Iuse)/ctx.L4CIp; rII=(ctx.C_II_mean*Iuse)/ctx.L4IIp;
    field=reshape(Euse,fieldRows,[]); filtered=conv2(padarray(field,[1 1],'circular'),ctx.L6Kernel,'same');
    l6=L6Convert(filtered(2:end-1,2:end-1),ctx.L6pars{ctx.L6parId}); l6=min(max(l6(:),3),3*numel(R.l6Mesh))./3;
    if epoch==targetEpoch, q=struct('recurrentE',[rES;rEC;rEI],'recurrentI',[rIS;rIC;rII],'l6',l6); return; end
    S0=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.S,rES,rIS,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    C0=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.C,rEC,rIC,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    I0=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.I,rEI,rII,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    S=ctx.p*S0+(1-ctx.p)*S; C=ctx.p*C0+(1-ctx.p)*C; I=ctx.p*I0+(1-ctx.p)*I;
end
end

function qr=make_queries(q,K,R)
n=numel(q.l6); classId=mod((0:n-1)',K)+1; phase=(1:K)'/max(K,1);
e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3);
baseE=median(e,1); baseI=median(i,1); baseL6=median(q.l6);
perturb=0.02*sin(2*pi*phase)+0.01*cos(6*pi*phase);
eClass=baseE.*(1+perturb); iClass=baseI.*(1-0.7*perturb);
l6Class=baseL6.*(1+0.5*perturb);
l6Class=min(max(l6Class,min(R.l6Mesh)),max(R.l6Mesh));
qr=struct('recurrentE',reshape(eClass(classId,:),[],1),'recurrentI', ...
    reshape(iClass(classId,:),[],1),'l6',l6Class(classId),'actualK',numel(unique(classId)));
end

function output=cg_pixel(q,ctx,cg)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3); output=zeros(n,3);
for p=1:3, output(:,p)=figure1g_cg_population_response(e(:,p),i(:,p),q.l6,ctx.PixLGNCtgr,cg.F(:,p)); end
end

function output=cg_class(q,ctx,cg)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3); output=zeros(n,3);
for p=1:3
    key=[e(:,p),i(:,p),q.l6]; [u,~,group]=unique(key,'rows'); response=zeros(size(u,1),5);
    for category=1:5, response(:,category)=cg.F{category,p}(u(:,2),u(:,1),u(:,3)); end
    output(:,p)=sum(ctx.PixLGNCtgr.*response(group,:),2);
end
end
