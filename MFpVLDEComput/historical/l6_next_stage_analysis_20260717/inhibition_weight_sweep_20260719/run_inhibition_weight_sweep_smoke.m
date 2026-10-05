function summary = run_inhibition_weight_sweep_smoke(sourceFile,outputRoot)
% Verify the inhibitory-weight pencil before full eigenspectrum jobs.

if nargin < 1 || isempty(sourceFile)
    sourceFile = getenv('IW_SOURCE_MAT');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = getenv('IW_OUT_ROOT');
end
if isempty(sourceFile) || ~isfile(sourceFile)
    error('L6NS:IWeightSource','IW_SOURCE_MAT must identify two_pathway_result.mat.');
end
if isempty(outputRoot)
    error('L6NS:IWeightOutput','IW_OUT_ROOT is required.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(sourceFile,'result');
pathway = loaded.result.Pathway;
jW0 = pathway.JRest+pathway.J6+pathway.JI;
reconstructionError = norm(jW0-pathway.JBaseline,'fro')/ ...
    max(norm(pathway.JBaseline,'fro'),eps);
assert(reconstructionError < 1e-14, ...
    'w_I=0 does not reconstruct the flexible-inhibition baseline.');

opts = struct('tol',1e-11,'maxit',1200,'p',80,'disp',0);
lambdaW0 = eigs(jW0,8,'largestreal',opts);
lambdaW1 = eigs(pathway.JRest+pathway.J6,8,'largestreal',opts);
maxRealW0 = max(real(lambdaW0));
maxRealW1 = max(real(lambdaW1));
expectedW0 = loaded.result.OppositeEffect.baselineLambdaReal;
expectedW1 = loaded.result.OppositeEffect.frozenILambdaReal;
assert(abs(maxRealW0-expectedW0) < 1e-8, ...
    'w_I=0 leading eigenvalue does not match the audited baseline.');
assert(abs(maxRealW1-expectedW1) < 1e-7, ...
    'w_I=1 leading eigenvalue does not match frozen inhibition.');

summary = table(0,1,reconstructionError,maxRealW0,expectedW0, ...
    maxRealW1,expectedW1,'VariableNames', ...
    {'fixedL6Weight','gamma6','baselineReconstructionError', ...
    'maxRealWI0','expectedMaxRealWI0','maxRealWI1','expectedMaxRealWI1'});
writetable(summary,fullfile(outputRoot,'inhibition_weight_smoke.tsv'), ...
    'FileType','text','Delimiter','\t');
fprintf('I-weight smoke passed: w_I=0 %.12g, w_I=1 %.12g.\n', ...
    maxRealW0,maxRealW1);
end
