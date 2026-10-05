function summary=run_all_spatial_controls_reequilibrated_task(taskIndex,setupFile,outputRoot)
% Re-equilibrate each spatial control before computing its Jacobian spectrum.

if nargin<1 || isempty(taskIndex)
    taskIndex=str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin<2 || isempty(setupFile); setupFile=getenv('SPATIAL_CONTROL_SETUP'); end
if nargin<3 || isempty(outputRoot); outputRoot=getenv('SPATIAL_CONTROL_OUTPUT'); end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

controls=["baseline";"shared_source_permutation"; ...
    "joint_alpha_0p10";"joint_alpha_0p25";"joint_alpha_0p50"; ...
    "remove_L6_smoothing";"remove_L4_smoothing"; ...
    "remove_L4_and_L6_smoothing";"L4_left_half";"L4_triangle"];
if taskIndex<1 || taskIndex>numel(controls)
    error('SpatialControl:Task','Task index must be between 1 and %d.',numel(controls));
end

loaded=load(setupFile,'setup');
setup=loaded.setup;
context=setup.Context;
control=controls(taskIndex);
[operators,metadata]=l6ns_full_spatial_control_operators(context,control);

initial=context.FixedPoint(:);
if control=="L4_left_half" || control=="L4_triangle"
    fprintf('%s: tracking the equilibrium from baseline by continuation.\n',control);
    [baselineOperators,~]=l6ns_full_spatial_control_operators(context,"baseline");
    [fixedPoint,continuationIterations,converged,continuationResidual, ...
        continuationSeconds]=local_continuation(initial,setup, ...
        baselineOperators,operators);
    iterations=continuationIterations;
    iterationSeconds=continuationSeconds;
    residualHistory=[iterations continuationResidual NaN];
    relaxation=NaN;
    solverMethod="continuation_fsolve_analytic_jacobian";
else
    [fixedPoint,residualHistory,iterations,converged,relaxation,solverMethod, ...
        iterationSeconds]=local_solve(initial,setup,operators);
end
fixedPointResidual=norm(l6ns_controlled_phi(fixedPoint,context,operators)-fixedPoint)/ ...
    max(norm(fixedPoint),eps);
relativeBaselineDistance=norm(fixedPoint-initial)/max(norm(initial),eps);
if ~converged
    warning('SpatialControl:NotConverged', ...
        '%s stopped at residual %.6e.',control,fixedPointResidual);
end

fprintf('%s: constructing analytic Jacobian at its own equilibrium.\n',control);
J=l6ns_controlled_jacobian(setup,operators,fixedPoint);
eigenvalues=eig(full(J),'vector');
singularValues=svd(full(J));
n=numel(fixedPoint)/3;
e=1:2*n; inhibitory=2*n+(1:n);
A=sparse(J(e,e)); B=sparse(J(e,inhibitory));
C=sparse(J(inhibitory,e)); D=sparse(J(inhibitory,inhibitory));
H=sparse(D-C*(A\B));
schurEigenvalues=eig(full(H),'vector');

baselineJacobianRelativeError=NaN;
if control=="baseline"
    baselineJacobianRelativeError=norm(J-setup.Pathway.JBaseline,'fro')/ ...
        max(norm(setup.Pathway.JBaseline,'fro'),eps);
end

summary=table(taskIndex,control,iterations,converged,relaxation,solverMethod, ...
    fixedPointResidual,iterationSeconds,relativeBaselineDistance, ...
    max(real(eigenvalues)),min(real(eigenvalues)), ...
    sum(abs(eigenvalues)<=0.01),sum(abs(eigenvalues)<=0.025), ...
    sum(abs(eigenvalues)<=0.05),sum(singularValues<=0.05), ...
    sum(abs(schurEigenvalues)<=0.05),baselineJacobianRelativeError, ...
    metadata.MaximumRowSumRelativeError,metadata.L6RowSumRelativeError, ...
    'VariableNames',{'taskIndex','control','iterations','converged', ...
    'relaxation','solverMethod','fixedPointResidual','iterationSeconds', ...
    'relativeBaselineDistance','maximumRealEigenvalue','minimumRealEigenvalue', ...
    'countAbsEigLE0p01','countAbsEigLE0p025','countAbsEigLE0p05', ...
    'countSingularLE0p05','countSchurEigLE0p05', ...
    'baselineJacobianRelativeError','maximumL4RowSumRelativeError', ...
    'l6RowSumRelativeError'});

tag=char(control);
writetable(summary,fullfile(outputRoot,[tag '.tsv']), ...
    'FileType','text','Delimiter','\t');
equilibrium=struct('S',fixedPoint(1:n),'C',fixedPoint(n+(1:n)), ...
    'I',fixedPoint(2*n+(1:n)));
save(fullfile(outputRoot,[tag '.mat']),'summary','equilibrium','fixedPoint', ...
    'residualHistory','eigenvalues','singularValues','schurEigenvalues', ...
    'metadata','-v7.3');
fprintf('%s: residual %.3e; eig/sigma/Schur <=.05: %d/%d/%d.\n',tag, ...
    fixedPointResidual,summary.countAbsEigLE0p05, ...
    summary.countSingularLE0p05,summary.countSchurEigLE0p05);
end

function [state,totalIterations,converged,residual,elapsed]= ...
        local_continuation(initial,setup,baselineOperators,targetOperators)
timer=tic;
state=initial;
totalIterations=0;
position=0;
step=0.10;
minimumStep=0.00625;
converged=false;
residual=Inf;
while position<1
    candidatePosition=min(1,position+step);
    operators=local_blend_operators(baselineOperators,targetOperators, ...
        candidatePosition);
    [candidate,exitflag,output,residual]=local_root( ...
        state,setup,operators,'none',80);
    totalIterations=totalIterations+output.iterations;
    fprintf('continuation %.5f residual %.6e exitflag %d\n', ...
        candidatePosition,residual,exitflag);
    if exitflag>0 && residual<=1e-9
        state=candidate;
        position=candidatePosition;
        step=min(0.15,step*1.25);
    else
        step=step/2;
        if step<minimumStep
            elapsed=toc(timer);
            return
        end
    end
end
converged=true;
elapsed=toc(timer);
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
    'FunctionTolerance',1e-10,'StepTolerance',1e-10, ...
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

function [state,residualHistory,iterations,converged,relaxation,solverMethod,elapsed]= ...
        local_solve(initial,setup,operators)
context=setup.Context;
tolerance=1e-9;
checkEvery=25;
relaxations=[0.33 0.15 0.07];
attemptLimits=[500 1000 2000];
elapsedTimer=tic;
bestState=initial;
bestResidual=Inf;
residualHistory=[];
iterations=0;
converged=false;
solverMethod="relaxed_iteration";

for attempt=1:numel(relaxations)
    relaxation=relaxations(attempt);
    state=bestState;
    previousCheck=Inf;
    for localIteration=1:attemptLimits(attempt)
        response=l6ns_controlled_phi(state,context,operators);
        if any(~isfinite(response))
            break
        end
        state=(1-relaxation)*state+relaxation*response;
        if mod(localIteration,checkEvery)==0 || localIteration==1
            residual=norm(l6ns_controlled_phi(state,context,operators)-state)/ ...
                max(norm(state),eps);
            residualHistory(end+1,:)=[iterations+localIteration residual relaxation]; %#ok<AGROW>
            fprintf('attempt %d iter %d p %.2f residual %.6e\n', ...
                attempt,iterations+localIteration,relaxation,residual);
            if residual<bestResidual
                bestResidual=residual;
                bestState=state;
            end
            if residual<=tolerance
                converged=true;
                iterations=iterations+localIteration;
                elapsed=toc(elapsedTimer);
                return
            end
            if residual>max(1e-3,previousCheck*20)
                break
            end
            previousCheck=residual;
        end
    end
    iterations=iterations+localIteration;
end

% An unstable fixed point cannot be reached by relaxed Picard iteration.
% Continue with a root solve of Phi(x)-x using its exact analytic Jacobian.
solverMethod="fsolve_analytic_jacobian";
relaxation=NaN;
[state,exitflag,output,residual]=local_root(bestState,setup,operators,'iter',150);
iterations=iterations+output.iterations;
residualHistory(end+1,:)=[iterations residual NaN];
converged=exitflag>0 && residual<=tolerance;
elapsed=toc(elapsedTimer);
end
