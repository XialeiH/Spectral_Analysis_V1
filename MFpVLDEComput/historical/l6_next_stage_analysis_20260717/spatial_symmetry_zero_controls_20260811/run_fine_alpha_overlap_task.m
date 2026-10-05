function summary = run_fine_alpha_overlap_task(taskIndex,setupFile,outputRoot)
% Resolve threshold crossings and spatial-smoothing/E-I overlap.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(setupFile)
    setupFile = getenv('SPATIAL_CONTROL_SETUP');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('SPATIAL_CONTROL_OUTPUT');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

alphaValues = [0 0.0025 0.005 0.01 0.02 0.03 0.05 0.075 0.10 ...
    0.125 0.15 0.20 0.25 0.35 0.50 0.75 1.00];
if taskIndex < 1 || taskIndex > numel(alphaValues)
    error('FineAlpha:Task','Task index must be 1--%d.',numel(alphaValues));
end
alpha = alphaValues(taskIndex);

loaded = load(setupFile,'setup');
setup = loaded.setup;
J = l6ns_spatial_smoothing_control_jacobian(setup,alpha,alpha);
n = size(J,1)/3;
e = 1:2*n;
i = 2*n+(1:n);
A = sparse(J(e,e));
B = sparse(J(e,i));
C = sparse(J(i,e));
D = sparse(J(i,i));

fprintf('Fine alpha %.4f: full spectrum and singular values.\n',alpha);
fullEigenvalues = eig(full(J),'vector');
fullSingularValues = svd(full(J));

% Removing B and C breaks the closed E-I feedback loop while retaining A and D.
fprintf('Fine alpha %.4f: decoupled E-I blocks.\n',alpha);
eigenvaluesA = eig(full(A),'vector');
eigenvaluesD = eig(full(D),'vector');
singularValuesA = svd(full(A));
singularValuesD = svd(full(D));
decoupledEigenvalues = [eigenvaluesA;eigenvaluesD];
decoupledSingularValues = [singularValuesA;singularValuesD];

% The low-lambda E-I branch is approximated by H v = lambda M v.
fprintf('Fine alpha %.4f: Schur generalized branch.\n',alpha);
X = A\B;
H = sparse(D-C*X);
M = sparse(speye(n)+C*(A\X));
schurEigenvalues = eig(full(H),'vector');
generalizedSchurEigenvalues = eig(full(H),full(M),'vector');

thresholds = [0.01 0.025 0.05];
fullEigCounts = local_counts(fullEigenvalues,thresholds);
fullSingularCounts = local_counts(fullSingularValues,thresholds);
decoupledEigCounts = local_counts(decoupledEigenvalues,thresholds);
decoupledSingularCounts = local_counts(decoupledSingularValues,thresholds);
schurCounts = local_counts(schurEigenvalues,thresholds);
generalizedSchurCounts = local_counts(generalizedSchurEigenvalues,thresholds);

summary = table(taskIndex,alpha, ...
    fullEigCounts(1),fullEigCounts(2),fullEigCounts(3), ...
    fullSingularCounts(1),fullSingularCounts(2),fullSingularCounts(3), ...
    decoupledEigCounts(1),decoupledEigCounts(2),decoupledEigCounts(3), ...
    decoupledSingularCounts(1),decoupledSingularCounts(2),decoupledSingularCounts(3), ...
    schurCounts(1),schurCounts(2),schurCounts(3), ...
    generalizedSchurCounts(1),generalizedSchurCounts(2),generalizedSchurCounts(3), ...
    norm(H,'fro')/(norm(D,'fro')+norm(C*X,'fro')),condest(A), ...
    'VariableNames',{'taskIndex','alpha', ...
    'fullEigLE0p01','fullEigLE0p025','fullEigLE0p05', ...
    'fullSigmaLE0p01','fullSigmaLE0p025','fullSigmaLE0p05', ...
    'decoupledEigLE0p01','decoupledEigLE0p025','decoupledEigLE0p05', ...
    'decoupledSigmaLE0p01','decoupledSigmaLE0p025','decoupledSigmaLE0p05', ...
    'schurEigLE0p01','schurEigLE0p025','schurEigLE0p05', ...
    'generalizedSchurEigLE0p01','generalizedSchurEigLE0p025', ...
    'generalizedSchurEigLE0p05','schurCancellationRatio','condestA'});

tag = strrep(sprintf('%.4f',alpha),'.','p');
writetable(summary,fullfile(outputRoot,['fine_alpha_' tag '.tsv']), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,['fine_alpha_' tag '.mat']), ...
    'summary','fullEigenvalues','fullSingularValues', ...
    'decoupledEigenvalues','decoupledSingularValues', ...
    'schurEigenvalues','generalizedSchurEigenvalues','-v7.3');
fprintf('alpha %.4f: full eig/sigma <=.05 %d/%d; Schur/generalized %d/%d.\n', ...
    alpha,fullEigCounts(3),fullSingularCounts(3), ...
    schurCounts(3),generalizedSchurCounts(3));
end

function counts = local_counts(values,thresholds)
counts = arrayfun(@(threshold) sum(abs(values)<=threshold),thresholds);
end
