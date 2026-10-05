function run_figure1g_batch_throughput(sourceRoot,outputRoot,fieldHC,repeatId)
% Measure recurrent and response throughput for simultaneous condition batches.
fieldHC=str2double(string(fieldHC)); repeatId=str2double(string(repeatId));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin'); addpath(fullfile(sourceRoot,'code'),'-begin'); addpath(fullfile(sourceRoot,'fast'),'-begin');
threads=str2double(getenv('SLURM_CPUS_PER_TASK')); if isfinite(threads)&&threads>=1,maxNumCompThreads(threads);end
B0=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
R=load(fullfile(sourceRoot,'bundles','cg_response_c100_a0.mat')); cg=figure1g_prepare_cg_interpolants(R);
base=h96_nn_prepare_base(B0.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat')); fast=h96_nn_set_condition_staticfirst(base,100,0);
n=fast.n; batchSizes=[1 4 16 64 256]; methods=["CG-Optimized-Pixel","DNN-Pixel"];
rows=cell(numel(batchSizes)*2,1); rowId=0;
for batchId=1:numel(batchSizes)
    batch=batchSizes(batchId); phase=reshape(0.01*sin(2*pi*(0:batch-1)/max(batch,1)),1,[]);
    S=B0.initialState.S.*(1+phase); C=B0.initialState.C.*(1-0.5*phase); I=B0.initialState.I.*(1+0.3*phase);
    E=.6923*S+.3077*C; scale=L6Convert(E,B0.ctx.EKpUse)./E; Suse=S.*scale; Cuse=C.*scale;
    Iuse=InhMulp(I,B0.ctx.IKpUse); recurrentFn=@() recurrent_batch(fast,Suse,Cuse,Iuse);
    [recurrentE,recurrentI]=recurrentFn(); l6=repmat(median(R.l6Mesh),n,batch);
    responseFns={@() cg_batch(recurrentE,recurrentI,l6,B0.ctx.PixLGNCtgr,cg), ...
        @() dnn_batch(recurrentE,recurrentI,l6,fast)};
    [recurrentMedian,~,~]=timing(@() recurrentFn());
    for methodId=1:2
        [responseMedian,responseMean,responseStd]=timing(responseFns{methodId}); total=recurrentMedian+responseMedian;
        rowId=rowId+1; rows{rowId}=table(fieldHC,repeatId,batch,methods(methodId),recurrentMedian,responseMedian,total, ...
            responseMean,responseStd,batch/max(total,eps),n*batch/max(total,eps), ...
            'VariableNames',{'FieldHC','Repeat','BatchSize','Method','RecurrentSeconds','ResponseSeconds', ...
            'CombinedSeconds','ResponseMeanSeconds','ResponseStdSeconds','FieldIterationsPerSecond','PixelUpdatesPerSecond'});
    end
end
T=vertcat(rows{:}); out=fullfile(outputRoot,'batch_throughput',sprintf('batch_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function [recurrentE,recurrentI]=recurrent_batch(fast,S,C,I)
recurrentE=fast.AE*[S;C]; recurrentI=fast.AI*I;
end

function output=cg_batch(recurrentE,recurrentI,l6,weights,cg)
n=size(l6,1); batch=size(l6,2); output=zeros(n,3,batch);
chunkSize=max(1,min(batch,floor(2e7/max(n*5,1))));
for first=1:chunkSize:batch
    ids=first:min(first+chunkSize-1,batch);
    for p=1:3
        e=recurrentE((p-1)*n+(1:n),ids); i=recurrentI((p-1)*n+(1:n),ids); response=zeros(n,numel(ids),5);
        for category=1:5, response(:,:,category)=reshape(cg.F{category,p}(i(:),e(:),l6(:,ids)),n,numel(ids)); end
        output(:,p,ids)=reshape(sum(response.*reshape(weights,n,1,5),3),n,1,numel(ids));
    end
end
end

function output=dnn_batch(recurrentE,recurrentI,l6Axis,fast)
n=fast.n; batch=size(l6Axis,2); rows=fast.mixRows; dynamicE=guard(reshape(recurrentE,n,3,batch)); dynamicI=reshape(recurrentI,n,3,batch);
pairCount=fast.activePairCount; output=zeros(n,3,batch); chunkSize=max(1,min(batch,floor(5e7/max(pairCount*96,1))));
for first=1:chunkSize:batch
    ids=first:min(first+chunkSize-1,batch); chunk=numel(ids);
    for p=1:3
        l6Rows=reshape(l6Axis(rows,ids),[],1); eRows=reshape(dynamicE(rows,p,ids),[],1); iRows=reshape(dynamicI(rows,p,ids),[],1);
        dynamic=zeros(pairCount*chunk,3,fast.nnPrecision);
        dynamic(:,1)=(cast(l6Rows,fast.nnPrecision)-fast.mu(1,3,p))./fast.sd(1,3,p);
        dynamic(:,2)=(cast(eRows,fast.nnPrecision)-fast.mu(1,4,p))./fast.sd(1,4,p);
        dynamic(:,3)=(cast(iRows,fast.nnPrecision)-fast.mu(1,5,p))./fast.sd(1,5,p);
        h=tanh(repmat(fast.staticLayer1(:,:,p),chunk,1)+dynamic*fast.dynamicW1T(:,:,p));
        h=tanh(h*fast.W2T(:,:,p)+fast.b2(:,:,p)); h=tanh(h*fast.W3T(:,:,p)+fast.b3(:,:,p));
        z=h*fast.W4T(:,:,p)+fast.b4(:,:,p); y=max(z,0)+log1p(exp(-abs(z)));
        y=reshape(y,pairCount,chunk); output(:,p,ids)=reshape(double(fast.mixAggregate*y),n,1,chunk);
    end
end
end

function x=guard(x)
edge=47500; full=55000; u=min(max((x-edge)./(full-edge),0),1); g=u.^3.*(10-15.*u+6.*u.^2); x=(1-g).*x+g.*(edge+0.02.*(x-edge));
end

function [med,avg,sd]=timing(fn)
fn(); samples=zeros(5,1); for k=1:numel(samples),timer=tic;fn();samples(k)=toc(timer);end; med=median(samples);avg=mean(samples);sd=std(samples);
end
