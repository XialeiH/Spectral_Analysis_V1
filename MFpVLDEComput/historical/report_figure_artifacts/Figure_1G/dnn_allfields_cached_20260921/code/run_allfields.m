function run_allfields(sourceRoot,outputRoot,fieldHC,repeatId,variant,phase)
% Cold iteration-only timing; all validation is performed after the timer.
fieldHC=str2double(string(fieldHC));
assert(ismember(fieldHC,[3 4 6 8 10 20 30 40]));
repeatId=str2double(string(repeatId)); variant=char(variant); phase=char(phase);
variants={'reference','paged_combined'};
assert(ismember(variant,variants));
addpath(fullfile(sourceRoot,'repo','Utils'),'-begin');
addpath(fullfile(sourceRoot,'code'),'-begin');
addpath(fullfile(sourceRoot,'fast'),'-begin');
addpath(fullfile(outputRoot,'code'),'-begin');
threads=8;
if strcmp(variant,'reference'), threads=16; end
maxNumCompThreads(threads);
assert(maxNumCompThreads==threads);
B=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
initial=B.initialState;
if repeatId>3
    rng(1729+repeatId,'twister');
    initial.S=max(0,initial.S.*exp(.08*randn(size(initial.S))-.5*.08^2));
    initial.C=max(0,initial.C.*exp(.08*randn(size(initial.C))-.5*.08^2));
    initial.I=max(0,initial.I.*exp(.08*randn(size(initial.I))-.5*.08^2));
end
prepTimer=tic;
base=h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
original=h96_nn_set_condition_symmetry(base,100,0);
fast=prepare_variant(original,variant);
prepareSeconds=toc(prepTimer);
step=@(x) update_state(x,fast,variant);
reference=@(x) update_state(x,original,'reference');
if strcmp(phase,'profile'), profile clear; profile on; end
[state,iterations,seconds,residual,converged]=iterate(step,initial);
if strcmp(phase,'profile')
    profile off;
    profileInfo=profile('info');
    save(fullfile(outputRoot,phase,['profile_' variant '.mat']),'profileInfo');
end

% Paired trajectories independently evolve from the same initial state.
candidate=initial; control=initial; stepMax=0; stepMean=0; stepResidualDifference=0;
for k=1:200
    a=step(candidate); b=reference(control);
    difference=abs([a.S-b.S;a.C-b.C;a.I-b.I]);
    stepMax=max(stepMax,max(difference)); stepMean=max(stepMean,mean(difference));
    stepResidualDifference=max(stepResidualDifference,abs(distance(a,candidate)-distance(b,control)));
    candidate=a; control=b;
    if k>=iterations, break; end
end
[referenceState,referenceIterations,~,referenceResidual,referenceConverged]=iterate(reference,initial);
delta=abs([state.S-referenceState.S;state.C-referenceState.C;state.I-referenceState.I]);
maxAbsolute=max(delta); meanAbsolute=mean(delta);
passed=converged && referenceConverged && iterations==referenceIterations && ...
    maxAbsolute<=1e-3 && meanAbsolute<=1e-4 && stepMax<=1e-3 && stepMean<=1e-4 && ...
    abs(residual-referenceResidual)<=1e-6 && stepResidualDifference<=1e-6;
allocatedCPUs=str2double(getenv('SLURM_CPUS_PER_TASK'));
node=string(getenv('SLURMD_NODENAME'));
T=table(fieldHC,repeatId,string(variant),seconds,iterations,residual,converged, ...
    threads,allocatedCPUs,node,prepareSeconds,numel(original.symRepresentative),fast.evalCount, ...
    referenceIterations,maxAbsolute,meanAbsolute,stepMax,stepMean, ...
    abs(residual-referenceResidual),stepResidualDifference,passed, ...
    'VariableNames',{'FieldHC','Repeat','Variant','Seconds','Iterations','FinalResidual', ...
    'Converged','NNThreads','AllocatedCPUs','Node','PrepareSeconds','OriginalRepresentatives', ...
    'EvaluatedInputs','ReferenceIterations','MaxAbsolute','MeanAbsolute','StepMaxAbsolute', ...
    'StepMeanAbsolute','ResidualDifference','StepResidualDifference','Passed'});
stem=sprintf('%s_%02dHC_repeat_%03d',variant,fieldHC,repeatId);
save(fullfile(outputRoot,phase,[stem '.mat']),'T','state','referenceState');
writetable(T,fullfile(outputRoot,phase,[stem '.tsv']),'FileType','text','Delimiter','\t');
fprintf('%s repeat=%d seconds=%.9f iterations=%d inputs=%d stepMax=%.9g passed=%d\n', ...
    variant,repeatId,seconds,iterations,fast.evalCount,stepMax,passed);
assert(passed,'Candidate failed step-by-step numerical validation.');
end

function fast=prepare_variant(original,variant)
fast=original;
rep=fast.symRepresentative;
pixels=fast.mixRows(rep);
static=fast.staticLayer1(rep,:,:);
inverse=(1:numel(rep)).';
if ~ismember(variant,{'reference','compact'})
    % Merge only identical pixel indices AND exactly identical static NN inputs.
    % This does not introduce any new spatial or numerical approximation.
    key=[double(pixels),double(reshape(static,numel(rep),[]))];
    [~,keep,inverse]=unique(key,'rows');
    pixels=pixels(keep); static=static(keep,:,:);
end
fast.repPixels=pixels;
fast.cachedStatic=static;
fast.outputGroup=inverse(fast.symGroup);
fast.evalCount=numel(pixels);
fast.selected=startsWith(variant,'selected_');
if fast.selected
    rows=pixels+(0:2)*fast.n;
    fast.AE=fast.AE(rows(:),:);
    fast.AI=fast.AI(rows(:),:);
end
fast.transposed=contains(variant,'transposed');
fast.paged=contains(variant,'paged');
fast.useQuad=ismember(variant,{'paged_quad','paged_combined'});
fast.usePadding=ismember(variant,{'paged_padding','paged_combined'});
fast.useAggregate=ismember(variant,{'paged_aggregate','paged_combined'});
if fast.useQuad
    assert(strcmp(fast.EKp{end},'quadratic'));
    xd=fast.EKp{4}; yd=fast.EKp{3};
    fast.quadCoeffs=[xd.^2,xd,ones(3,1)]\yd;
end
if fast.usePadding
    fast.padRows=[fast.fieldRows,1:fast.fieldRows,1];
    fast.padCols=[fast.fieldCols,1:fast.fieldCols,1];
end
if fast.useAggregate
    [ii,jj,vv]=find(fast.mixAggregate);
    fast.collapsedAggregate=cast(sparse(ii,fast.outputGroup(jj),double(vv), ...
        fast.n,fast.evalCount),'like',fast.mixAggregate);
end
if fast.transposed
    fast.cachedStatic=permute(fast.cachedStatic,[2 1 3]);
    fast.dynamicW1T=permute(fast.dynamicW1T,[2 1 3]);
    fast.W2T=permute(fast.W2T,[2 1 3]);
    fast.W3T=permute(fast.W3T,[2 1 3]);
    fast.W4T=permute(fast.W4T,[2 1 3]);
    fast.b2=permute(fast.b2,[2 1 3]);
    fast.b3=permute(fast.b3,[2 1 3]);
end
end

function [state,iterations,seconds,residual,converged]=iterate(step,initial)
state=initial; stableCount=0; converged=false;
timer=tic;
for iterations=1:200
    next=step(state);
    residual=distance(next,state);
    state=next;
    if residual<5e-3, stableCount=stableCount+1; else, stableCount=0; end
    if stableCount>=3, converged=true; break; end
end
seconds=toc(timer);
end

function next=update_state(state,fast,variant)
S=state.S; C=state.C; I=state.I;
E=.6923*S+.3077*C;
if isfield(fast,'useQuad') && fast.useQuad
    pars=fast.EKp; coeff=fast.quadCoeffs;
    converted=pars{1}*(E-pars{2});
    xAdj=E(E>pars{3}(1));
    a=coeff(1); b=coeff(2); c=coeff(3);
    converted(converted>pars{4}(1))=(-b+sqrt(b^2-4*a*c+4*a*xAdj))/(2*a);
    scale=converted./E;
else
    scale=L6Convert(E,fast.EKp)./E;
end
Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp);
Euse=.6923*Suse+.3077*Cuse;
recurrentE=fast.AE*[Suse;Cuse]; recurrentI=fast.AI*Iuse;
field=reshape(Euse,fast.fieldRows,fast.fieldCols);
if isfield(fast,'usePadding') && fast.usePadding
    padded=field(fast.padRows,fast.padCols);
else
    padded=padarray(field,[1 1],'circular');
end
filtered=conv2(padded,fast.kernel,'same');
l6Input=filtered(2:end-1,2:end-1);
if isempty(fast.l6PP), l6=L6Convert(l6Input,fast.l6Pars); else, l6=ppval(fast.l6PP,l6Input); end
l6=min(max(l6(:),3),fast.l6Max)./3;
if strcmp(variant,'reference')
    output=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fast);
else
    output=response_compact(recurrentE,recurrentI,l6,fast);
end
next=struct('S',fast.p*output(:,1)+fast.oneMinusP*S, ...
    'C',fast.p*output(:,2)+fast.oneMinusP*C, ...
    'I',fast.p*output(:,3)+fast.oneMinusP*I);
end

function output=response_compact(recurrentE,recurrentI,l6,fast)
count=fast.evalCount; pixels=fast.repPixels;
if fast.selected
    dynamicE=guard(reshape(recurrentE,count,3));
    dynamicI=reshape(recurrentI,count,3);
else
    e=reshape(recurrentE,fast.n,3); i=reshape(recurrentI,fast.n,3);
    dynamicE=guard(e(pixels,:)); dynamicI=i(pixels,:);
end
dynamic=zeros(count,3,3,'single');
l6Rows=single(l6(pixels));
for population=1:3
    dynamic(:,1,population)=(l6Rows-fast.mu(1,3,population))./fast.sd(1,3,population);
    dynamic(:,2,population)=(single(dynamicE(:,population))-fast.mu(1,4,population))./fast.sd(1,4,population);
    dynamic(:,3,population)=(single(dynamicI(:,population))-fast.mu(1,5,population))./fast.sd(1,5,population);
end
if fast.transposed, dynamic=permute(dynamic,[2 1 3]); end
if fast.paged
    if fast.transposed
        h=tanh(fast.cachedStatic+pagemtimes(fast.dynamicW1T,dynamic));
        h=tanh(pagemtimes(fast.W2T,h)+fast.b2);
        h=tanh(pagemtimes(fast.W3T,h)+fast.b3);
        z=pagemtimes(fast.W4T,h)+fast.b4;
        z=permute(z,[2 1 3]);
    else
        h=tanh(fast.cachedStatic+pagemtimes(dynamic,fast.dynamicW1T));
        h=tanh(pagemtimes(h,fast.W2T)+fast.b2);
        h=tanh(pagemtimes(h,fast.W3T)+fast.b3);
        z=pagemtimes(h,fast.W4T)+fast.b4;
    end
else
    z=zeros(count,1,3,'single');
    for population=1:3
        if fast.transposed
            h=tanh(fast.cachedStatic(:,:,population)+fast.dynamicW1T(:,:,population)*dynamic(:,:,population));
            h=tanh(fast.W2T(:,:,population)*h+fast.b2(:,:,population));
            h=tanh(fast.W3T(:,:,population)*h+fast.b3(:,:,population));
            zp=fast.W4T(:,:,population)*h+fast.b4(:,:,population);
            z(:,:,population)=zp.';
        else
            h=tanh(fast.cachedStatic(:,:,population)+dynamic(:,:,population)*fast.dynamicW1T(:,:,population));
            h=tanh(h*fast.W2T(:,:,population)+fast.b2(:,:,population));
            h=tanh(h*fast.W3T(:,:,population)+fast.b3(:,:,population));
            z(:,:,population)=h*fast.W4T(:,:,population)+fast.b4(:,:,population);
        end
    end
end
y=max(z,0)+log1p(exp(-abs(z)));
output=zeros(fast.n,3);
for population=1:3
    if fast.useAggregate
        output(:,population)=double(fast.collapsedAggregate*y(:,1,population));
    else
        output(:,population)=double(fast.mixAggregate*y(fast.outputGroup,1,population));
    end
end
end

function x=guard(x)
edge=47500; full=55000;
u=min(max((x-edge)./(full-edge),0),1);
g=u.^3.*(10-15.*u+6.*u.^2);
x=(1-g).*x+g.*(edge+.02.*(x-edge));
end

function d=distance(a,b)
d=norm([a.S-b.S;a.C-b.C;a.I-b.I])/max(norm([b.S;b.C;b.I]),eps);
end
