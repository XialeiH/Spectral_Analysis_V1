function summary = run_inhibition_weight_eigenspectrum_task(taskIndex,sourceFile,outputRoot)
% Compute one full eigenspectrum for fixed dynamic L6 and one w_I value.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceFile)
    sourceFile = getenv('IW_SOURCE_MAT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('IW_OUT_ROOT');
end
values = l6ns_inhibition_weight_values();
if ~isscalar(taskIndex) || ~isfinite(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(values)
    error('L6NS:IWeightTask','Task index must be an integer from 1 to %d.',numel(values));
end
if isempty(sourceFile) || ~isfile(sourceFile)
    error('L6NS:IWeightSource','IW_SOURCE_MAT must identify two_pathway_result.mat.');
end
if isempty(outputRoot)
    error('L6NS:IWeightOutput','IW_OUT_ROOT is required.');
end
dataRoot = fullfile(outputRoot,'eigdata');
if ~exist(dataRoot,'dir'); mkdir(dataRoot); end

wI = values(taskIndex);
gammaI = 1-wI;
tag = l6ns_inhibition_weight_tag(wI);
outputFile = fullfile(dataRoot,sprintf('inhibition_weight_eigvals_wI%s.mat',tag));
summaryFile = fullfile(dataRoot,sprintf('inhibition_weight_summary_wI%s.tsv',tag));

loaded = load(sourceFile,'result');
pathway = loaded.result.Pathway;
jacobian = pathway.JRest+pathway.J6+gammaI*pathway.JI;
baselineReconstructionError = NaN;
if wI == 0
    baselineReconstructionError = norm(jacobian-pathway.JBaseline,'fro')/ ...
        max(norm(pathway.JBaseline,'fro'),eps);
    assert(baselineReconstructionError < 1e-14, ...
        'w_I=0 does not reconstruct the baseline Jacobian.');
end

fprintf('Computing full eigenspectrum for w_I=%+.2f, gamma_I=%+.2f.\n',wI,gammaI);
timer = tic;
eigenvalues = eig(full(jacobian));
runtimeSeconds = toc(timer);
assert(numel(eigenvalues)==size(jacobian,1) && all(isfinite(eigenvalues)), ...
    'Full eigenspectrum is incomplete or nonfinite.');

realPart = real(eigenvalues);
sortedReal = sort(realPart);
summary = table(wI,gammaI,numel(eigenvalues),min(realPart), ...
    local_percentile(sortedReal,1),local_percentile(sortedReal,5), ...
    median(realPart),local_percentile(sortedReal,95), ...
    local_percentile(sortedReal,99),max(realPart),max(abs(eigenvalues)), ...
    max(abs(imag(eigenvalues))),max(realPart)<1,runtimeSeconds, ...
    baselineReconstructionError,string(outputFile), ...
    'VariableNames',{'wI','gammaI','numberEigenvalues','minReal','p01Real', ...
    'p05Real','medianReal','p95Real','p99Real','maxReal','spectralRadius', ...
    'maxAbsImag','stable','runtimeSeconds','baselineReconstructionError','matFile'});
metadata = struct('Definition','J=Jrest+J6+(1-wI)*JI', ...
    'FixedL6Weight',0,'Gamma6',1,'WI',wI,'GammaI',gammaI, ...
    'SourceFile',sourceFile,'TaskIndex',taskIndex);
save(outputFile,'eigenvalues','summary','metadata','-v7.3');
writetable(summary,summaryFile,'FileType','text','Delimiter','\t');
fprintf('Saved %d eigenvalues for w_I=%+.2f in %.1f seconds.\n', ...
    numel(eigenvalues),wI,runtimeSeconds);
end

function value = local_percentile(sortedValues,percent)
position = 1+(numel(sortedValues)-1)*percent/100;
lower = floor(position);
upper = ceil(position);
fraction = position-lower;
value = (1-fraction)*sortedValues(lower)+fraction*sortedValues(upper);
end
