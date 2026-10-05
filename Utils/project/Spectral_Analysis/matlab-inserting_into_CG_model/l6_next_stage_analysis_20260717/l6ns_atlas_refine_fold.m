function event = l6ns_atlas_refine_fold(setup,pathwayName,stateA,wA,stateB,wB,outputDir,options)
% Refine a real unit-eigenvalue event along a pseudo-arclength segment.

if nargin<8; options=struct(); end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
targetTolerance=local_option(options,'TargetTolerance',1e-10);
maximumIterations=local_option(options,'MaximumIterations',36);
etaW=local_option(options,'EtaW',10);

fixed=setup.Context.FixedPoint(:);
n=numel(fixed)/3;
population=reshape(fixed,n,3);
populationScales=max(sqrt(mean(population.^2,1)),1);
metricWeights=repelem(1./(n*populationScales.^2),n)';
phi=local_phi(setup,pathwayName);
[tangentState,tangentW,totalDistance]=local_secant(stateA,wA,stateB,wB, ...
    metricWeights,etaW);

alphaA=local_alpha(setup,pathwayName,stateA,wA);
alphaB=local_alpha(setup,pathwayName,stateB,wB);
if (alphaA-1)*(alphaB-1)>0
    error('Input states do not bracket max Re lambda = 1.');
end
lowFraction=0;
highFraction=1;
lowState=stateA(:); lowW=wA; lowAlpha=alphaA;
highState=stateB(:); highW=wB; highAlpha=alphaB;
iterationRows={};

for iteration=1:maximumIterations
    fraction=(lowFraction+highFraction)/2;
    guessState=(1-fraction)*stateA+fraction*stateB;
    guessW=(1-fraction)*wA+fraction*wB;
    [state,w,corrector]=local_corrector(setup,pathwayName,phi,stateA,wA, ...
        tangentState,tangentW,fraction*totalDistance,guessState,guessW, ...
        metricWeights,etaW,fixed);
    if ~corrector.converged
        error('Arclength event corrector failed at fraction %.12g.',fraction);
    end
    alpha=local_alpha(setup,pathwayName,state,w);
    iterationRows(end+1,:)={iteration,fraction,w,alpha,alpha-1, ... %#ok<AGROW>
        corrector.residualNorm,corrector.relativeResidual};
    if abs(alpha-1)<targetTolerance
        lowState=state; lowW=w; lowAlpha=alpha;
        highState=state; highW=w; highAlpha=alpha;
        break
    end
    if (lowAlpha-1)*(alpha-1)<=0
        highFraction=fraction; highState=state; highW=w; highAlpha=alpha;
    else
        lowFraction=fraction; lowState=state; lowW=w; lowAlpha=alpha;
    end
end

if abs(lowAlpha-1)<=abs(highAlpha-1)
    state=lowState; w=lowW; alpha=lowAlpha;
else
    state=highState; w=highW; alpha=highAlpha;
end
J=l6ns_state_jacobian(setup,pathwayName,state,w);
modes=l6ns_eigenpairs(J,24,'largestreal',setup.Config);
[~,criticalIndex]=min(abs(modes.Lambda-1));
lambda=modes.Lambda(criticalIndex);
r=real(modes.Right(:,criticalIndex));
r=r/norm(r);
l=real(modes.Left(:,criticalIndex));
l=l/(l'*r);

parameterStep=2e-6*max(1,abs(w));
Rw=(phi(state,w+parameterStep)-phi(state,w-parameterStep))/(2*parameterStep);
etaSN=real(l'*Rw);
stateScale=max(1,norm(state)/sqrt(numel(state)));
steps=stateScale*[1e-3 3e-4 1e-4];
bValues=nan(size(steps));
center=phi(state,w)-state;
for index=1:numel(steps)
    h=steps(index);
    plus=phi(state+h*r,w)-(state+h*r);
    minus=phi(state-h*r,w)-(state-h*r);
    Brr=(plus-2*center+minus)/(h^2);
    bValues(index)=0.5*real(l'*Brr);
end
bSN=median(bValues);
other=modes.Lambda;
other(criticalIndex)=[];
spectralSeparation=min(abs(other-lambda));
kernelDimension=sum(abs(modes.Lambda-1)<1e-6);
diagnostic=l6ns_atlas_point_metrics(setup,pathwayName,state,w,tangentW, ...
    struct('ModeCount',40,'MaximumModeCount',80));
classification="unresolved_real_singularity";
if abs(imag(lambda))<1e-9 && kernelDimension==1 && ...
        abs(etaSN)>1e-8 && abs(bSN)>1e-8
    classification="generic_saddle_node";
end
summary=table(string(pathwayName),w,alpha,real(lambda),imag(lambda), ...
    etaSN,bSN,spectralSeparation,kernelDimension,min(abs(bValues)),max(abs(bValues)), ...
    string(classification), ...
    'VariableNames',{'pathway','freezeWeight','maxRealLambda','lambdaReal', ...
    'lambdaImag','etaSN','bSN','spectralSeparation','kernelDimension', ...
    'minimumAbsBEstimate','maximumAbsBEstimate','classification'});
iterations=cell2table(iterationRows,'VariableNames',{'iteration','segmentFraction', ...
    'freezeWeight','maxRealLambda','boundaryResidual','residualNorm','relativeResidual'});
coefficientTable=table(steps(:),bValues(:),'VariableNames',{'step','bSN'});
event=struct('Summary',summary,'State',state,'Right',r,'Left',l, ...
    'PointDiagnostic',diagnostic,'Iterations',iterations,'CoefficientConvergence',coefficientTable);
writetable(summary,fullfile(outputDir,'inhibition_fold_classification.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(iterations,fullfile(outputDir,'inhibition_fold_refinement.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(coefficientTable,fullfile(outputDir,'inhibition_fold_coefficient_convergence.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,'inhibition_fold_event.mat'),'event','-v7.3');
end

function alpha=local_alpha(setup,pathwayName,state,w)
J=l6ns_state_jacobian(setup,pathwayName,state,w);
modes=l6ns_eigenpairs(J,8,'largestreal',setup.Config);
alpha=max(real(modes.Lambda));
end

function [tangentState,tangentW,distance]=local_secant(stateA,wA,stateB,wB,metricWeights,etaW)
tangentState=stateB(:)-stateA(:);
tangentW=wB-wA;
distance=sqrt((metricWeights.*tangentState)'*tangentState+etaW^2*tangentW^2);
tangentState=tangentState/distance;
tangentW=tangentW/distance;
end

function [state,w,diagnostic]=local_corrector(setup,pathwayName,phi,baseState,baseW, ...
        tangentState,tangentW,targetDistance,state,w,metricWeights,etaW,fixed)
tolerance=2e-10*max(1,norm(fixed));
converged=false;
residualNorm=inf;
for iteration=1:18
    residualState=phi(state,w)-state;
    arcResidual=(metricWeights.*tangentState)'*(state-baseState)+ ...
        etaW^2*tangentW*(w-baseW)-targetDistance;
    residual=[residualState;arcResidual];
    residualNorm=norm(residual);
    if residualNorm<tolerance; converged=true; break; end
    A=l6ns_state_jacobian(setup,pathwayName,state,w)-speye(numel(state));
    parameterStep=2e-6*max(1,abs(w));
    Rw=(phi(state,w+parameterStep)-phi(state,w-parameterStep))/(2*parameterStep);
    bordered=[A,sparse(Rw);sparse((metricWeights.*tangentState)'),etaW^2*tangentW];
    delta=bordered\(-residual);
    if any(~isfinite(delta)); break; end
    accepted=false;
    lineScale=1;
    for lineIteration=1:10
        candidateState=state+lineScale*delta(1:end-1);
        candidateW=w+lineScale*delta(end);
        candidateResidual=phi(candidateState,candidateW)-candidateState;
        candidateArc=(metricWeights.*tangentState)'*(candidateState-baseState)+ ...
            etaW^2*tangentW*(candidateW-baseW)-targetDistance;
        if norm([candidateResidual;candidateArc])<residualNorm
            state=candidateState; w=candidateW; accepted=true; break
        end
        lineScale=lineScale/2;
    end
    if ~accepted; break; end
end
diagnostic=struct('converged',converged,'iterations',iteration, ...
    'residualNorm',residualNorm,'relativeResidual',residualNorm/max(1,norm(state)));
end

function phi=local_phi(setup,pathwayName)
if strcmpi(pathwayName,'L6')
    phi=@(x,w)l6ns_phi(x,w,setup.Context,[1 1],1);
else
    phi=@(x,w)l6ns_phi(x,0,setup.Context,[1 1],1-w);
end
end

function value=local_option(options,name,defaultValue)
if isfield(options,name); value=options.(name); else; value=defaultValue; end
end
