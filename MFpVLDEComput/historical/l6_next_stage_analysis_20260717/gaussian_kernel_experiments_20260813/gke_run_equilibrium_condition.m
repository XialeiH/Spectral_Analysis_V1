function summary=gke_run_equilibrium_condition(taskIndex,setupFile,outputRoot)
% Solve and save one manipulation-specific equilibrium before spectral work.

if nargin<1 || isempty(taskIndex); taskIndex=str2double(getenv('SLURM_ARRAY_TASK_ID')); end
if nargin<2 || isempty(setupFile); setupFile=getenv('GKE_SETUP'); end
if nargin<3 || isempty(outputRoot); outputRoot=getenv('GKE_OUTPUT'); end
catalog=gke_condition_catalog();
condition=catalog(taskIndex,:);
dataRoot=fullfile(outputRoot,'equilibrium_only','data');
summaryRoot=fullfile(outputRoot,'equilibrium_only','summaries');
if ~exist(dataRoot,'dir'); mkdir(dataRoot); end
if ~exist(summaryRoot,'dir'); mkdir(summaryRoot); end

loaded=load(setupFile,'setup');
setup=loaded.setup;
context=setup.Context;
[operators,metadata,operatorAudit]=gke_build_operators(context,condition.name);
knownFoldControls=["swap_EI_ranges" "remove_L6_smoothing"];
[fixedPoint,solver]=gke_track_equilibrium(setup,operators, ...
    condition.name=="baseline",ismember(condition.name,knownFoldControls));
response=l6ns_controlled_phi(fixedPoint,context,operators);
fixedPointResidual=norm(response-fixedPoint)/max(norm(fixedPoint),eps);
if ~solver.Converged || fixedPointResidual>1e-9
    error('GKE:FixedPoint','%s did not reach a valid fixed point.',condition.name);
end

n=numel(fixedPoint)/3;
equilibrium=struct('S',fixedPoint(1:n),'C',fixedPoint(n+(1:n)), ...
    'I',fixedPoint(2*n+(1:n)));
equilibrium.EWeighted=(1-context.CWeight)*equilibrium.S+ ...
    context.CWeight*equilibrium.C;
equilibrium.ESum=equilibrium.S+equilibrium.C;
summary=table(condition.id,condition.name,condition.group,condition.label, ...
    fixedPointResidual,solver.Iterations,solver.Seconds,min(fixedPoint), ...
    max(fixedPoint),norm(fixedPoint-context.FixedPoint(:))/ ...
    max(norm(context.FixedPoint(:)),eps), ...
    'VariableNames',{'id','name','group','label','fixedPointResidual', ...
    'solverIterations','solverSeconds','minimumFiringRate', ...
    'maximumFiringRate','relativeBaselineDistance'});
save(fullfile(dataRoot,condition.name+'.mat'),'condition','metadata', ...
    'operatorAudit','equilibrium','fixedPoint','solver','summary','-v7.3');
writetable(summary,fullfile(summaryRoot,condition.name+'.tsv'), ...
    'FileType','text','Delimiter','\t');
fprintf('%s equilibrium complete: residual %.3e, range [%.6f, %.6f].\n', ...
    condition.name,fixedPointResidual,min(fixedPoint),max(fixedPoint));
end
