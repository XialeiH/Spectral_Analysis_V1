function run_figure1g_response_replay(sourceRoot,outputRoot,fieldHC,repeatId)
% Replay realistic response queries through CG and h96 evaluators only.
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
interpolants = figure1g_prepare_cg_interpolants(R);

snapshots = record_cg_queries(B.initialState,B.ctx,R,[1 10 25 50]);
methods = ["CG-Paper3-Pixel","CG-Paper3-Class", ...
    "CG-Optimized-Pixel","CG-Optimized-Class","DNN-Pixel","DNN-Class"];
rows = cell(numel(snapshots)*numel(methods),1);
rowId = 0;
for snapshotId = 1:numel(snapshots)
    q = snapshots(snapshotId);
    evaluators = { ...
        @() cg_paper3_pixel(q,B.ctx,R), ...
        @() cg_paper3_class(q,B.ctx,R), ...
        @() cg_optimized_pixel(q,B.ctx,interpolants), ...
        @() cg_optimized_class(q,B.ctx,interpolants), ...
        @() h96_nn_all_population_responses_staticfirst(q.recurrentE,q.recurrentI,q.l6,fast), ...
        @() h96_nn_all_population_responses_unique(q.recurrentE,q.recurrentI,q.l6,fast)};
    reference = evaluators{1}();
    for methodId = 1:numel(methods)
        fn = evaluators{methodId};
        value = fn(); % warm-up outside timing
        samples = zeros(5,1);
        for timingId = 1:numel(samples)
            timer = tic; fn(); samples(timingId) = toc(timer);
        end
        rowId = rowId+1;
        rows{rowId} = table(fieldHC,repeatId,q.epoch,methods(methodId), ...
            median(samples),mean(samples),std(samples),q.queryCount,q.uniqueQueryCount, ...
            norm(value-reference,'fro')/max(norm(reference,'fro'),eps), ...
            'VariableNames',{'FieldHC','Repeat','Epoch','Method','MedianSeconds', ...
            'MeanSeconds','StdSeconds','QueryCount','UniqueQueryCount','RelativeToPaper3'});
    end
end
T = vertcat(rows{:});
out = fullfile(outputRoot,'response_replay', ...
    sprintf('replay_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function snapshots = record_cg_queries(initialState,ctx,R,targetEpochs)
S=initialState.S; C=initialState.C; I=initialState.I;
fieldRows=ctx.N_HCOutY*ctx.NPixY; pars=ctx.L6pars{ctx.L6parId};
snapshots = repmat(struct('epoch',0,'recurrentE',[],'recurrentI',[], ...
    'l6',[],'queryCount',0,'uniqueQueryCount',0),numel(targetEpochs),1);
snapshotId = 0;
for epoch=1:max(targetEpochs)
    E=.6923*S+.3077*C; adjusted=L6Convert(E,ctx.EKpUse); scale=adjusted./E;
    Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,ctx.IKpUse); Euse=.6923*Suse+.3077*Cuse;
    rES=(ctx.C_SS_meanU*Suse+ctx.C_SC_meanU*Cuse)/ctx.L4SEp;
    rEC=(ctx.C_CS_meanU*Suse+ctx.C_CC_meanU*Cuse)/ctx.L4CEp;
    rEI=(ctx.C_IS_mean*Suse+ctx.C_IC_mean*Cuse)/ctx.L4IEp;
    rIS=(ctx.C_SI_mean*Iuse)/ctx.L4SIp;
    rIC=(ctx.C_CI_mean*Iuse)/ctx.L4CIp;
    rII=(ctx.C_II_mean*Iuse)/ctx.L4IIp;
    field=reshape(Euse,fieldRows,[]); padded=padarray(field,[1 1],'circular');
    filtered=conv2(padded,ctx.L6Kernel,'same'); l6=L6Convert(filtered(2:end-1,2:end-1),pars);
    l6=min(max(l6(:),3),3*numel(R.l6Mesh))./3;
    if any(epoch==targetEpochs)
        snapshotId=snapshotId+1;
        recurrentE=[rES;rEC;rEI]; recurrentI=[rIS;rIC;rII];
        key=[recurrentE,recurrentI,repmat(l6,3,1)];
        snapshots(snapshotId)=struct('epoch',epoch,'recurrentE',recurrentE, ...
            'recurrentI',recurrentI,'l6',l6,'queryCount',size(key,1), ...
            'uniqueQueryCount',size(unique(key,'rows'),1));
    end
    Sout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.S,rES,rIS,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    Cout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.C,rEC,rIC,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    Iout=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.I,rEI,rII,ctx.PixLGNCtgr,l6,[1;2;3;4]);
    S=ctx.p*Sout+(1-ctx.p)*S; C=ctx.p*Cout+(1-ctx.p)*C; I=ctx.p*Iout+(1-ctx.p)*I;
end
end

function output = cg_paper3_pixel(q,ctx,R)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3);
output=zeros(n,3);
output(:,1)=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.S,e(:,1),i(:,1),ctx.PixLGNCtgr,q.l6,[1;2;3;4]);
output(:,2)=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.C,e(:,2),i(:,2),ctx.PixLGNCtgr,q.l6,[1;2;3;4]);
output(:,3)=LDEIterFunc_Grating_135Func_RealLGNL6(R.l4EMesh,R.l4IMesh,R.l6Mesh,R.func.I,e(:,3),i(:,3),ctx.PixLGNCtgr,q.l6,[1;2;3;4]);
end

function output = cg_paper3_class(q,ctx,R)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3);
output=zeros(n,3); funcs={R.func.S,R.func.C,R.func.I};
for population=1:3
    key=[e(:,population),i(:,population),q.l6];
    [uniqueKey,~,group]=unique(key,'rows');
    weights=zeros(size(uniqueKey,1),5);
    [XX,YY,ZZ]=meshgrid(unique(R.l4EMesh),unique(R.l4IMesh),R.l6Mesh);
    for category=1:5
        weights(:,category)=interp3(XX,YY,ZZ,squeeze(funcs{population}(category,:,:,:)), ...
            uniqueKey(:,1),uniqueKey(:,2),uniqueKey(:,3),'linear');
    end
    output(:,population)=sum(ctx.PixLGNCtgr.*weights(group,:),2);
end
end

function output = cg_optimized_pixel(q,ctx,interpolants)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3);
output=zeros(n,3);
for population=1:3
    output(:,population)=figure1g_cg_population_response(e(:,population),i(:,population), ...
        q.l6,ctx.PixLGNCtgr,interpolants.F(:,population));
end
end

function output = cg_optimized_class(q,ctx,interpolants)
n=numel(q.l6); e=reshape(q.recurrentE,n,3); i=reshape(q.recurrentI,n,3);
output=zeros(n,3);
for population=1:3
    key=[e(:,population),i(:,population),q.l6];
    [uniqueKey,~,group]=unique(key,'rows');
    responses=zeros(size(uniqueKey,1),5);
    for category=1:5
        responses(:,category)=interpolants.F{category,population}(uniqueKey(:,2),uniqueKey(:,1),uniqueKey(:,3));
    end
    output(:,population)=sum(ctx.PixLGNCtgr.*responses(group,:),2);
end
end
