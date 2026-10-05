function event = l6ns_atlas_refine_sector_event(setup,pathwayName,stateA,wA, ...
        stateB,wB,rowWave,columnWave,outputStem,options)
% Localize a real unit-eigenvalue event in one translation sector.

if nargin<10; options=struct(); end
targetTolerance=local_option(options,'TargetTolerance',1e-9);
maximumIterations=local_option(options,'MaximumIterations',30);
etaW=local_option(options,'EtaW',10);
fixed=setup.Context.FixedPoint(:);
n=numel(fixed)/3;
population=reshape(fixed,n,3);
populationScales=max(sqrt(mean(population.^2,1)),1);
metricWeights=repelem(1./(n*populationScales.^2),n)';
phi=local_phi(setup,pathwayName);
Q=l6ns_atlas_sector_basis(setup.Context.MapSize,rowWave,columnWave);
[tangentState,tangentW,totalDistance]=local_secant(stateA,wA,stateB,wB, ...
    metricWeights,etaW);
alphaA=local_sector_alpha(setup,pathwayName,stateA,wA,Q);
alphaB=local_sector_alpha(setup,pathwayName,stateB,wB,Q);
if (alphaA-1)*(alphaB-1)>0
    error('Input states do not bracket a unit eigenvalue in sector (%d,%d).', ...
        rowWave,columnWave);
end

lowFraction=0; highFraction=1;
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
        error('Sector event corrector failed at fraction %.12g.',fraction);
    end
    alpha=local_sector_alpha(setup,pathwayName,state,w,Q);
    iterationRows(end+1,:)={iteration,fraction,w,alpha,alpha-1, ... %#ok<AGROW>
        corrector.residualNorm,corrector.relativeResidual};
    if abs(alpha-1)<targetTolerance
        lowState=state; lowW=w; lowAlpha=alpha;
        highState=state; highW=w; highAlpha=alpha;
        break
    elseif (lowAlpha-1)*(alpha-1)<=0
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
reduced=Q'*J*Q;
modes=l6ns_eigenpairs(reduced,16,'largestreal',setup.Config);
[~,criticalIndex]=min(abs(modes.Lambda-1));
lambda=modes.Lambda(criticalIndex);
r=Q*modes.Right(:,criticalIndex);
l=Q*modes.Left(:,criticalIndex);
phase=exp(-1i*angle(r'*real(r)));
r=real(r*phase); r=r/norm(r);
l=real(l*phase); l=l/(l'*r);
hParameter=2e-6*max(1,abs(w));
Rw=(phi(state,w+hParameter)-phi(state,w-hParameter))/(2*hParameter);
etaSN=real(l'*Rw);
stateScale=max(1,norm(state)/sqrt(numel(state)));
steps=stateScale*[1e-3 5e-4 2e-4 1e-4];
bValues=nan(size(steps));
center=phi(state,w)-state;
for index=1:numel(steps)
    h=steps(index);
    plus=phi(state+h*r,w)-(state+h*r);
    minus=phi(state-h*r,w)-(state-h*r);
    Brr=(plus-2*center+minus)/(h^2);
    bValues(index)=0.5*real(l'*Brr);
end
bSN=median(bValues(end-2:end));
other=modes.Lambda; other(criticalIndex)=[];
spectralSeparation=min(abs(other-lambda));
classification="unresolved_sector_event";
if abs(imag(lambda))<1e-8 && abs(etaSN)>1e-8 && abs(bSN)>1e-8
    classification="generic_saddle_node_in_sector";
end
summary=table(string(pathwayName),w,rowWave,columnWave,alpha,real(lambda), ...
    imag(lambda),etaSN,bSN,spectralSeparation,string(classification), ...
    'VariableNames',{'pathway','freezeWeight','rowWaveNumber', ...
    'columnWaveNumber','sectorMaximumRealEigenvalue','lambdaReal','lambdaImag', ...
    'etaSN','bSN','sectorSpectralSeparation','classification'});
iterations=cell2table(iterationRows,'VariableNames',{'iteration','segmentFraction', ...
    'freezeWeight','sectorMaximumRealEigenvalue','boundaryResidual', ...
    'residualNorm','relativeResidual'});
coefficientTable=table(steps(:),bValues(:),'VariableNames',{'step','bSN'});
event=struct('Summary',summary,'State',state,'Right',r,'Left',l, ...
    'Iterations',iterations,'CoefficientConvergence',coefficientTable);
[outputDir,stem]=fileparts(outputStem);
if ~exist(outputDir,'dir'); mkdir(outputDir); end
save(fullfile(outputDir,[stem '.mat']),'event','-v7.3');
writetable(summary,fullfile(outputDir,[stem '_summary.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(iterations,fullfile(outputDir,[stem '_refinement.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(coefficientTable,fullfile(outputDir,[stem '_coefficients.tsv']), ...
    'FileType','text','Delimiter','\t');
end

function alpha=local_sector_alpha(setup,pathwayName,state,w,Q)
J=l6ns_state_jacobian(setup,pathwayName,state,w);
lambda=eigs(Q'*J*Q,4,'largestreal');
alpha=max(real(lambda));
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
converged=false; residualNorm=inf;
for iteration=1:18
    residualState=phi(state,w)-state;
    arcResidual=(metricWeights.*tangentState)'*(state-baseState)+ ...
        etaW^2*tangentW*(w-baseW)-targetDistance;
    residual=[residualState;arcResidual];
    residualNorm=norm(residual);
    if residualNorm<tolerance; converged=true; break; end
    A=l6ns_state_jacobian(setup,pathwayName,state,w)-speye(numel(state));
    h=2e-6*max(1,abs(w));
    Rw=(phi(state,w+h)-phi(state,w-h))/(2*h);
    bordered=[A,sparse(Rw);sparse((metricWeights.*tangentState)'),etaW^2*tangentW];
    delta=bordered\(-residual);
    if any(~isfinite(delta)); break; end
    accepted=false; lineScale=1;
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
