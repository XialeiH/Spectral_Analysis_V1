function result = l6ns_atlas_weighted_continue(setup,pathwayName,seedState,seedW,direction,outputFile,options)
% Bidirectional-ready weighted pseudo-arclength continuation from one root.

if nargin<7; options=struct(); end
minimumStep=local_option(options,'MinimumStep',2e-4);
maximumStep=local_option(options,'MaximumStep',0.08);
step=local_option(options,'InitialStep',0.01);
maximumPoints=local_option(options,'MaximumPoints',120);
weightBounds=local_option(options,'WeightBounds',[-0.5 1.2]);
etaW=local_option(options,'EtaW',10);
modeCount=local_option(options,'ModeCount',40);

fixed=setup.Context.FixedPoint(:);
n=numel(fixed)/3;
fixedPopulation=reshape(fixed,n,3);
populationScales=max(sqrt(mean(fixedPopulation.^2,1)),1);
metricWeights=repelem(1./(n*populationScales.^2),n)';
phi=local_phi(setup,pathwayName);

[seedState,seedDiagnostic]=l6ns_fixed_point_newton(phi,fixed,seedW,seedState(:));
if ~seedDiagnostic.converged
    error('Seed polishing failed at w=%.12g.',seedW);
end
[tangentState,tangentW]=local_regular_tangent(setup,pathwayName,seedState,seedW, ...
    direction,metricWeights,etaW);

states={seedState};
summaryRows={};
modeTables={};
modeVectors={};
seedMetrics=l6ns_atlas_point_metrics(setup,pathwayName,seedState,seedW,tangentW, ...
    struct('ModeCount',modeCount,'MaximumModeCount',max(2*modeCount,80)));
summaryRows{1}=seedMetrics.Summary;
modeTables{1}=seedMetrics.Modes;
modeVectors{1}=struct('Right',seedMetrics.RightVectors,'Left',seedMetrics.LeftVectors);
previousState=seedState;
previousW=seedW;
terminalReason="maximum_points";
attemptRows={};

for pointIndex=1:maximumPoints
    retries=0;
    accepted=false;
    while ~accepted && step>=minimumStep && retries<8
        guessState=previousState+step*tangentState;
        guessW=previousW+step*tangentW;
        [state,w,corrector]=local_corrector(setup,pathwayName,phi,previousState, ...
            previousW,tangentState,tangentW,step,guessState,guessW, ...
            metricWeights,etaW,fixed);
        retries=retries+1;
        attemptRows(end+1,:)={pointIndex,retries,step,w,corrector.converged, ... %#ok<AGROW>
            corrector.iterations,corrector.residualNorm,corrector.relativeResidual};
        if corrector.converged
            accepted=true;
        else
            step=step/2;
        end
    end
    if ~accepted
        terminalReason="unresolved_after_step_reduction";
        break
    end

    [newTangentState,newTangentW]=local_secant_tangent(previousState,previousW, ...
        state,w,tangentState,tangentW,metricWeights,etaW);
    states{end+1,1}=state; %#ok<AGROW>
    metrics=l6ns_atlas_point_metrics(setup,pathwayName,state,w,newTangentW, ...
        struct('ModeCount',modeCount,'MaximumModeCount',max(2*modeCount,80)));
    summaryRows{end+1,1}=metrics.Summary; %#ok<AGROW>
    modeTables{end+1,1}=metrics.Modes; %#ok<AGROW>
    modeVectors{end+1,1}=struct('Right',metrics.RightVectors, ... %#ok<AGROW>
        'Left',metrics.LeftVectors);

    previousState=state;
    previousW=w;
    tangentState=newTangentState;
    tangentW=newTangentW;
    if corrector.iterations<=4
        step=min(maximumStep,1.25*step);
    elseif corrector.iterations>=9
        step=max(minimumStep,0.65*step);
    end

    if w<weightBounds(1) || w>weightBounds(2)
        terminalReason="parameter_boundary";
        break
    elseif any(~isfinite(state)) || max(abs(state))>5000
        terminalReason="numerical_divergence_boundary";
        break
    elseif pointIndex>12 && ...
            norm(state-seedState)/max(1,norm(seedState))<1e-5 && ...
            abs(w-seedW)<1e-4
        terminalReason="closed_full_state_loop";
        break
    end
end

summary=vertcat(summaryRows{:});
attempts=cell2table(attemptRows,'VariableNames',{'pointIndex','retryIndex', ...
    'attemptedStep','freezeWeight','converged','newtonIterations', ...
    'residualNorm','relativeResidual'});
result=struct('Pathway',string(pathwayName),'Direction',direction, ...
    'Summary',summary,'States',{states},'ModeTables',{modeTables}, ...
    'ModeVectors',{modeVectors},'Attempts',attempts, ...
    'PopulationScales',populationScales,'EtaW',etaW,'TerminalReason',terminalReason);
if ~isempty(outputFile)
    outputDir=fileparts(outputFile);
    if ~exist(outputDir,'dir'); mkdir(outputDir); end
    save(outputFile,'result','-v7.3');
    [~,stem]=fileparts(outputFile);
    writetable(summary,fullfile(outputDir,[stem '_summary.tsv']), ...
        'FileType','text','Delimiter','\t');
    writetable(attempts,fullfile(outputDir,[stem '_attempts.tsv']), ...
        'FileType','text','Delimiter','\t');
end
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
    if residualNorm<tolerance
        converged=true;
        break
    end
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
        candidateResidualState=phi(candidateState,candidateW)-candidateState;
        candidateArc=(metricWeights.*tangentState)'*(candidateState-baseState)+ ...
            etaW^2*tangentW*(candidateW-baseW)-targetDistance;
        if norm([candidateResidualState;candidateArc])<residualNorm
            state=candidateState;
            w=candidateW;
            accepted=true;
            break
        end
        lineScale=lineScale/2;
    end
    if ~accepted; break; end
end
diagnostic=struct('converged',converged,'iterations',iteration, ...
    'residualNorm',residualNorm,'relativeResidual',residualNorm/max(1,norm(state)));
end

function [tangentState,tangentW]=local_regular_tangent(setup,pathwayName,state,w, ...
        direction,metricWeights,etaW)
phi=local_phi(setup,pathwayName);
A=l6ns_state_jacobian(setup,pathwayName,state,w)-speye(numel(state));
h=2e-6*max(1,abs(w));
Rw=(phi(state,w+h)-phi(state,w-h))/(2*h);
tangentW=sign(direction);
tangentState=A\(-Rw*tangentW);
[tangentState,tangentW]=local_normalize(tangentState,tangentW,metricWeights,etaW);
end

function [tangentState,tangentW]=local_secant_tangent(firstState,firstW,secondState, ...
        secondW,oldState,oldW,metricWeights,etaW)
tangentState=secondState-firstState;
tangentW=secondW-firstW;
[tangentState,tangentW]=local_normalize(tangentState,tangentW,metricWeights,etaW);
alignment=(metricWeights.*oldState)'*tangentState+etaW^2*oldW*tangentW;
if alignment<0; tangentState=-tangentState; tangentW=-tangentW; end
end

function [statePart,weightPart]=local_normalize(statePart,weightPart,metricWeights,etaW)
lengthValue=sqrt((metricWeights.*statePart)'*statePart+etaW^2*weightPart^2);
statePart=statePart/lengthValue;
weightPart=weightPart/lengthValue;
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
