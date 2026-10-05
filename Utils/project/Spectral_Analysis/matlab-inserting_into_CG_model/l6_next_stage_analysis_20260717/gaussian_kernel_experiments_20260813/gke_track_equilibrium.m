function [fixedPoint,solver]=gke_track_equilibrium( ...
        setup,targetOperators,isBaseline,startOnTargetBranch)
% Track the baseline-connected fixed-point branch with exact D(Phi-I).

if nargin<4; startOnTargetBranch=false; end

context=setup.Context;
[baselineOperators,~,~]=gke_build_operators(context,"baseline");
initial=context.FixedPoint(:);
timer=tic;
if isBaseline
    [fixedPoint,exitflag,output,residual]=local_root( ...
        initial,setup,baselineOperators,'iter',80);
    continuationPosition=1;
    continuationResidual=residual;
    continuationExitflag=exitflag;
    totalIterations=output.iterations;
else
    [state,exitflag,output,residual]=local_root( ...
        initial,setup,baselineOperators,'none',50);
    if exitflag<=0 || residual>1e-9
        error('GKE:BaselineRoot','Could not refine the baseline fixed point.');
    end
    totalIterations=output.iterations;
    if startOnTargetBranch
        [fixedPoint,exitflag,fallbackOutput,residual,fallbackMethod]= ...
            local_target_fallback(initial,initial,setup,targetOperators);
        totalIterations=totalIterations+fallbackOutput.iterations;
        continuationPosition=[];
        continuationResidual=[];
        continuationExitflag=[];
        solver=local_solver_struct(fixedPoint,residual,totalIterations,timer, ...
            continuationPosition,continuationResidual,continuationExitflag, ...
            "target branch requested after previously observed continuation fold; "+ ...
            fallbackMethod);
        return
    end
    position=0;
    step=0.10;
    minimumStep=0.003125;
    continuationPosition=[];
    continuationResidual=[];
    continuationExitflag=[];
    continuationFailed=false;
    while position<1
        candidatePosition=min(1,position+step);
        operators=local_blend_operators(baselineOperators,targetOperators,candidatePosition);
        [candidate,candidateExitflag,candidateOutput,candidateResidual]= ...
            local_root(state,setup,operators,'none',100);
        totalIterations=totalIterations+candidateOutput.iterations;
        continuationPosition(end+1,1)=candidatePosition; %#ok<AGROW>
        continuationResidual(end+1,1)=candidateResidual; %#ok<AGROW>
        continuationExitflag(end+1,1)=candidateExitflag; %#ok<AGROW>
        fprintf('continuation %.5f residual %.6e exitflag %d\n', ...
            candidatePosition,candidateResidual,candidateExitflag);
        if candidateExitflag>0 && candidateResidual<=1e-9 && min(candidate)>-1e-8
            state=candidate;
            position=candidatePosition;
            step=min(0.15,step*1.25);
        else
            step=step/2;
            if step<minimumStep
                continuationFailed=true;
                fprintf(['Baseline-connected continuation stopped at %.5f; ' ...
                    'trying target-system iteration seeds.\n'],position);
                break
            end
        end
    end
    if continuationFailed
        [fixedPoint,exitflag,fallbackOutput,residual,fallbackMethod]= ...
            local_target_fallback(initial,state,setup,targetOperators);
        totalIterations=totalIterations+fallbackOutput.iterations;
    else
        fixedPoint=state;
        residual=continuationResidual(end);
        fallbackMethod="none";
    end
end

method="analytic-Jacobian continuation from baseline";
if ~isBaseline && fallbackMethod~="none"
    method=method+"; "+fallbackMethod;
end
solver=local_solver_struct(fixedPoint,residual,totalIterations,timer, ...
    continuationPosition,continuationResidual,continuationExitflag,method);
end

function solver=local_solver_struct(~,residual,totalIterations,timer, ...
        continuationPosition,continuationResidual,continuationExitflag,method)
solver=struct();
solver.Method=method;
solver.Converged=residual<=1e-9;
solver.FixedPointResidual=residual;
solver.Iterations=totalIterations;
solver.Seconds=toc(timer);
solver.ContinuationPosition=continuationPosition;
solver.ContinuationResidual=continuationResidual;
solver.ContinuationExitflag=continuationExitflag;
end

function [state,exitflag,output,residual,method]= ...
        local_target_fallback(baselineSeed,branchSeed,setup,targetOperators)
% Find a target-system equilibrium when the baseline-connected branch folds.

context=setup.Context;
seedNames=strings(0,1);
seeds=cell(0,1);
relaxations=[0.5 0.2 0.1 1];
for seedIndex=1:2
    if seedIndex==1; base=baselineSeed; baseName="baseline";
    else; base=branchSeed; baseName="last branch"; end
    for relaxation=relaxations
        candidate=base;
        for iteration=1:500
            response=l6ns_controlled_phi(candidate,context,targetOperators);
            candidate=(1-relaxation)*candidate+relaxation*response;
            if any(~isfinite(candidate)) || norm(candidate)>1e8
                break
            end
            if mod(iteration,10)==0
                candidateResidual=norm( ...
                    l6ns_controlled_phi(candidate,context,targetOperators)-candidate)/ ...
                    max(norm(candidate),eps);
                if candidateResidual<=1e-8; break; end
            end
        end
        candidateResidual=norm( ...
            l6ns_controlled_phi(candidate,context,targetOperators)-candidate)/ ...
            max(norm(candidate),eps);
        fprintf('target iteration %-43s residual %.6e\n', ...
            sprintf('%s relaxation %.2f',baseName,relaxation),candidateResidual);
        if candidateResidual<=1e-9 && min(candidate)>-1e-8
            state=candidate;
            exitflag=1;
            output=struct('iterations',iteration);
            residual=candidateResidual;
            method=sprintf('target branch from %s, relaxation %.2f', ...
                baseName,relaxation);
            return
        end
        seeds{end+1,1}=candidate; %#ok<AGROW>
        seedNames(end+1,1)=sprintf('%s, 500-step relaxation %.2f', ...
            baseName,relaxation); %#ok<AGROW>
    end
end
seeds=[seeds;{baselineSeed};{branchSeed}];
seedNames=[seedNames;"baseline direct";"last branch direct"];

bestResidual=Inf;
bestState=branchSeed;
bestExitflag=-1;
bestOutput=struct('iterations',0);
for index=1:numel(seeds)
    seed=seeds{index};
    if any(~isfinite(seed)) || norm(seed)>1e8; continue; end
    [candidate,candidateExitflag,candidateOutput,candidateResidual]= ...
        local_root(seed,setup,targetOperators,'none',300);
    fprintf('target fallback %-42s residual %.6e exitflag %d\n', ...
        seedNames(index),candidateResidual,candidateExitflag);
    if candidateResidual<bestResidual && min(candidate)>-1e-8
        bestResidual=candidateResidual;
        bestState=candidate;
        bestExitflag=candidateExitflag;
        bestOutput=candidateOutput;
    end
    if candidateExitflag>0 && candidateResidual<=1e-9 && min(candidate)>-1e-8
        state=candidate;
        exitflag=candidateExitflag;
        output=candidateOutput;
        residual=candidateResidual;
        method="target branch from "+seedNames(index);
        return
    end
end
error('GKE:TargetRootFailed', ...
    'No target equilibrium found; best admissible residual %.6e (exitflag %d).', ...
    bestResidual,bestExitflag);
end

function operators=local_blend_operators(baseline,target,weight)
names={'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II','L6'};
for index=1:numel(names)
    name=names{index};
    operators.(name)=(1-weight)*baseline.(name)+weight*target.(name);
end
end

function [state,exitflag,output,residual]=local_root(initial,setup,operators,display,maxIterations)
context=setup.Context;
identity=speye(numel(initial));
options=optimoptions('fsolve','Algorithm','trust-region', ...
    'SpecifyObjectiveGradient',true,'Display',display, ...
    'FunctionTolerance',1e-11,'StepTolerance',1e-11, ...
    'MaxIterations',maxIterations,'MaxFunctionEvaluations',2*maxIterations);
[state,~,exitflag,output]=fsolve(@objective,initial,options);
residual=norm(l6ns_controlled_phi(state,context,operators)-state)/ ...
    max(norm(state),eps);

    function [value,jacobian]=objective(candidate)
        value=l6ns_controlled_phi(candidate,context,operators)-candidate;
        if nargout>1
            jacobian=l6ns_controlled_jacobian(setup,operators,candidate)-identity;
        end
    end
end
