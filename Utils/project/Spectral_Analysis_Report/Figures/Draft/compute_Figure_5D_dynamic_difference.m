function compute_Figure_5D_dynamic_difference(taskIds)
% Paired, fixed-input linear response comparison; no singular-mode optimization.
projectRoot = [repro_paths('project') ''];
modelRoot = fullfile(projectRoot,'Spectral_Analysis','matlab-inserting_into_CG_model');
analysisRoot = fullfile(modelRoot,'l6_next_stage_analysis_20260717');
artifactRoot = fullfile(modelRoot,'report_figure_artifacts', ...
    '6_Pathway_Compensation','Figure_5D');
outputRoot = fullfile(artifactRoot,'fixed_input_dynamic_difference');
if ~isfolder(outputRoot), mkdir(outputRoot); end
repro_addpath(fullfile(modelRoot,'Complete_Code_for_Paper3','NYU-Vision-2Drive-main','Utils'));
repro_addpath(analysisRoot,'-begin');
repro_addpath(fullfile(analysisRoot,'real_tuning_l6_i_20260722'),'-begin');
runtimeRoot = fullfile(analysisRoot,'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722','runtime_h96');
repro_addpath(runtimeRoot,'-begin');
assert(startsWith(which('predict_pref6D_S'),runtimeRoot));
assert(startsWith(which('LocalResponse_6D_MLP_prefAngle'),runtimeRoot));
data = sortrows(readtable(fullfile(artifactRoot,'validation_s0025', ...
    'valley_metrics_s0025.tsv'),'FileType','text','Delimiter','\t'),'taskId');
if nargin < 1, taskIds = 1:height(data); end
loaded = load(fullfile(analysisRoot,'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat'),'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
n = numel(baseline)/3;
assert(n == 1600 && all(context.OrientationUse == 0,'all') && ...
    all(context.ContrastUse == 100,'all'));
tauMs = 10.3402405296839;
timesMs = 0:0.5:180;
initialHC = 1;
sourceFile = fullfile(modelRoot,'report_figure_artifacts','4_Pathway_Modes', ...
    'Figure_4H','data','figure_04_trajectory.mat');
source = load(sourceFile,'states');
initialDelta = double(source.states(:,1))-baseline;
clear source
initialDelta = initialDelta * (initialHC/sqrt(hc_squared(initialDelta,n)));
assert(abs(sqrt(hc_squared(initialDelta,n))-initialHC)<1e-12);
sourceDirectionFile = sourceFile;
J0 = real_tuning_true_jacobian(baseline,context,1,1);
baselineResidual = norm(l6ns_phi(baseline,0,context,[1 1],1,'true')-baseline)/norm(baseline);
assert(baselineResidual<1e-9,'Baseline no longer matches the preserved runtime.');
[baselineDelta,baselineSolverError] = response(J0,initialDelta,tauMs,timesMs,n,true);
baselineEnergy = trapz(timesMs,hc_squared(baselineDelta,n));
baselineHalfGridEnergy = trapz(timesMs(1:2:end),hc_squared(baselineDelta(:,1:2:end),n));
assert(baselineEnergy>0);
save(fullfile(outputRoot,'baseline.mat'),'initialDelta','initialHC','timesMs', ...
    'tauMs','baselineDelta','baselineEnergy','baselineSolverError','baselineResidual', ...
    'sourceDirectionFile');
fprintf('BASELINE residual %.3g; solver discrepancy %.3g; J nnz %d\n', ...
    baselineResidual,baselineSolverError,nnz(J0));
for taskId = taskIds
    timer = tic;
    beta6 = data.beta6(taskId); betaI = data.betaI(taskId);
    phi = @(state) l6ns_phi(state,-beta6,context,[1 1],1+betaI,'true');
    fixed = real_tuning_fixed_point(phi,baseline,context.RelaxationP);
    assert(fixed.Converged,'Fixed point did not converge.');
    staticError = sqrt(hc_squared(fixed.State-baseline,n));
    assert(abs(staticError-data.staticResponseError(taskId))<1e-5, ...
        'Static response does not match the original Figure 5D at task %d.',taskId);
    if taskId == 1
        delta = baselineDelta; solverError = baselineSolverError;
        derivativeError = 0;
    else
        J = real_tuning_true_jacobian(fixed.State,context,1+beta6,1+betaI);
        direction = initialDelta/norm(initialDelta);
        step = 1e-3;
        finiteDifference = (phi(fixed.State+step*direction)- ...
            phi(fixed.State-step*direction))/(2*step);
        derivativeError = norm(finiteDifference-J*direction)/max(norm(J*direction),eps);
        assert(derivativeError<1e-4,'Jacobian directional derivative check failed.');
        [delta,solverError] = response(J,initialDelta,tauMs,timesMs,n, ...
            ismember(taskId,[5 16 height(data)]));
    end
    differenceSquared = hc_squared(delta-baselineDelta,n);
    responseSquared = hc_squared(delta,n);
    dynamicDifference = sqrt(trapz(timesMs,differenceSquared)/baselineEnergy);
    halfGridDifference = sqrt(trapz(timesMs(1:2:end), ...
        differenceSquared(1:2:end))/baselineHalfGridEnergy);
    quadratureError = abs(dynamicDifference-halfGridDifference);
    assert(quadratureError<1e-3,'Time quadrature has not converged.');
    responseGain = sqrt(trapz(timesMs,responseSquared)/(timesMs(end)*initialHC^2));
    row = table(taskId,beta6,betaI,data.arcLength(taskId),staticError, ...
        dynamicDifference,responseGain,quadratureError,solverError,derivativeError, ...
        fixed.Residual,'VariableNames',{'taskId','beta6','betaI','arcLength', ...
        'staticResponseError','dynamicDifference','responseGain','quadratureError', ...
        'solverError','derivativeError','fixedPointResidual'});
    writetable(row,fullfile(outputRoot,sprintf('point_%02d.tsv',taskId)), ...
        'FileType','text','Delimiter','\t');
    fixedPoint = fixed.State;
    save(fullfile(outputRoot,sprintf('point_%02d.mat',taskId)), ...
        'row','fixedPoint','delta','timesMs','initialDelta','initialHC','tauMs');
    fprintf('TASK %02d static %.8g D_dyn %.8g G_response %.8g; %.1f s\n', ...
        taskId,staticError,dynamicDifference,responseGain,toc(timer));
end
files = arrayfun(@(k) fullfile(outputRoot,sprintf('point_%02d.tsv',k)), ...
    1:height(data),'UniformOutput',false);
if all(cellfun(@isfile,files))
    rows = cellfun(@(f) readtable(f,'FileType','text','Delimiter','\t'), ...
        files,'UniformOutput',false);
    results = vertcat(rows{:});
    writetable(results,fullfile(outputRoot,'dynamic_difference.tsv'), ...
        'FileType','text','Delimiter','\t');
end
end

function squared = hc_squared(delta,n)
e = 0.6923*delta(1:n,:)+0.3077*delta(n+(1:n),:);
i = delta(2*n+(1:n),:);
squared = mean(0.8*e.^2+0.2*i.^2,1);
end

function [delta,solverError] = response(J,initialDelta,tauMs,timesMs,n,verify)
A = (J-speye(size(J,1)))/tauMs;
options = odeset('RelTol',1e-8,'AbsTol',1e-11,'MaxStep',0.5);
solution = ode45(@(~,x) A*x,[timesMs(1) timesMs(end)],initialDelta,options);
delta = deval(solution,timesMs);
solverError = NaN;
if verify
    options = odeset('RelTol',1e-10,'AbsTol',1e-13,'MaxStep',0.25);
    refined = ode45(@(~,x) A*x,[timesMs(1) timesMs(end)],initialDelta,options);
    deltaRefined = deval(refined,timesMs);
    solverError = sqrt(trapz(timesMs,hc_squared(delta-deltaRefined,n))/ ...
        trapz(timesMs,hc_squared(deltaRefined,n)));
    assert(solverError<1e-5,'Linear ODE integration has not converged.');
end
end
