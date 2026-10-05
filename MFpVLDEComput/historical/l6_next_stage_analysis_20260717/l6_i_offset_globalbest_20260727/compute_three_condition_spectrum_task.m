function compute_three_condition_spectrum_task()
% Compute full spectra and leading eigen/singular modes for one comparison case.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
setupFile = getenv('THREE_SETUP_FILE');
orangeFile = getenv('THREE_ORANGE_FILE');
blueFile = getenv('THREE_BLUE_FILE');
outputRoot = getenv('THREE_OUTPUT_ROOT');
if ~isfinite(taskId) || taskId<1 || taskId>3 || ...
        ~isfile(setupFile) || ~isfile(orangeFile) || ~isfile(blueFile) || ...
        isempty(outputRoot)
    error('ThreeSpectrum:Environment','Invalid task or missing input paths.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
beta6 = [0;0.034;0.034];
betaI = [0;0.0460053091161357;0.0176002942418064];
caseLabel = ["baseline"; ...
    "orange FPP coordinate, true moved state"; ...
    "blue HC-minimum coordinate, true moved state"];

switch taskId
    case 1
        state = context.FixedPoint(:);
        J = setup.Pathway.JBaseline;
    case 2
        one = load(orangeFile,'fixed');
        state = one.fixed.State(:);
        J = real_tuning_true_jacobian(state,context,1+beta6(taskId),1+betaI(taskId));
    case 3
        one = load(blueFile,'fixedPoint');
        state = one.fixedPoint(:);
        J = real_tuning_true_jacobian(state,context,1+beta6(taskId),1+betaI(taskId));
end

timer = tic;
eigenvalues = eig(full(J),'vector');
fullSpectrumSeconds = toc(timer);
[rightMode,leadingLambda] = local_leading_right(J);
leftMode = local_matching_left(J,leadingLambda,rightMode);
svdsOptions = struct('tol',1e-9,'maxit',2500,'p',80,'disp',0);
[singularOutput,singularValue,singularInput] = svds(J,1,'largest',svdsOptions);
rightMode = local_real_phase(rightMode);
leftMode = local_real_phase(leftMode);
singularInput = local_real_phase(singularInput);
singularOutput = local_real_phase(singularOutput);

summary = table(taskId,caseLabel(taskId),beta6(taskId),betaI(taskId), ...
    real(leadingLambda),imag(leadingLambda),max(real(eigenvalues)), ...
    min(real(eigenvalues)),max(abs(imag(eigenvalues))),singularValue, ...
    fullSpectrumSeconds,numel(eigenvalues), ...
    'VariableNames',{'taskId','caseLabel','beta6','betaI', ...
    'leadingReal','leadingImag','maximumReal','minimumReal', ...
    'maximumAbsoluteImaginary','largestSingularValue', ...
    'fullSpectrumSeconds','numberEigenvalues'});
outputFile = fullfile(outputRoot,sprintf('condition_%d_spectrum_modes.mat',taskId));
save(outputFile,'summary','state','eigenvalues','rightMode','leftMode', ...
    'singularInput','singularOutput','singularValue','-v7.3');
writetable(summary,fullfile(outputRoot,sprintf('condition_%d_summary.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
disp(summary);
end

function [mode,lambda] = local_leading_right(J)
options = struct('tol',1e-10,'maxit',2500,'p',100,'isreal',true,'disp',0);
try
    [vectors,values] = eigs(J,10,'largestreal',options);
catch
    [vectors,values] = eigs(J,10,'lr',options);
end
lambdaPool = diag(values);
[~,index] = max(real(lambdaPool));
lambda = lambdaPool(index);
mode = vectors(:,index)/norm(vectors(:,index));
end

function mode = local_matching_left(J,lambda,rightMode)
options = struct('tol',1e-10,'maxit',2500,'p',100,'isreal',true,'disp',0);
try
    [vectors,values] = eigs(J',10,'largestreal',options);
catch
    [vectors,values] = eigs(J',10,'lr',options);
end
lambdaPool = diag(values);
[~,index] = min(abs(lambdaPool-conj(lambda)));
mode = vectors(:,index);
overlap = mode'*rightMode;
if abs(overlap)>1e-12
    mode = mode/conj(overlap);
end
mode = mode/norm(mode);
end

function vector = local_real_phase(vector)
[~,index] = max(abs(vector));
vector = vector*exp(-1i*angle(vector(index)));
vector = real(vector);
vector = vector/max(norm(vector),eps);
end
