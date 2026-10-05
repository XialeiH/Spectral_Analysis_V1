function slopes = audit_fpp_theoretical_slopes(setupFile, outputRoot)
% Audit FPP top-eigenvalue sensitivities for L6, excitation, and inhibition.

if nargin < 1 || isempty(setupFile)
    codeRoot = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');
    setupFile = fullfile(codeRoot, 'results_global_bifurcation_20260722', ...
        'global_bifurcation_setup.mat');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile( ...
        [repro_paths('project') '/Spectral_Analysis'], ...
        'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
        'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors', ...
        'L6 and Inhibition');
end

loaded = load(setupFile, 'setup');
pathway = loaded.setup.Pathway;
j6 = pathway.J6;
jE = pathway.JRest;
jI = pathway.JI;
jBaseline = pathway.JBaseline;
n = pathway.PopulationSize;

[lambda, left, right] = local_leading_mode(jBaseline, loaded.setup.Config);
s6 = real(left' * (j6 * right));
sE = real(left' * (jE * right));
sI = real(left' * (jI * right));

step = 1e-5;
fd6 = local_central_sensitivity(jBaseline, j6, step, loaded.setup.Config);
fdE = local_central_sensitivity(jBaseline, jE, step, loaded.setup.Config);
fdI = local_central_sensitivity(jBaseline, jI, step, loaded.setup.Config);

slopes = struct();
slopes.L6Inhibition = -s6 / sI;
slopes.L6Excitation = -s6 / sE;
slopes.ExcitationInhibition = -sE / sI;
slopes.SensitivityL6 = s6;
slopes.SensitivityExcitation = sE;
slopes.SensitivityInhibition = sI;

reconstructionError = norm(jE + j6 + jI - jBaseline, 'fro') / ...
    max(norm(jBaseline, 'fro'), eps);
l6ISourceFraction = norm(j6(:, 2*n+(1:n)), 'fro') / max(norm(j6, 'fro'), eps);
excitationISourceFraction = norm(jE(:, 2*n+(1:n)), 'fro') / ...
    max(norm(jE, 'fro'), eps);
inhibitionNonISourceFraction = norm(jI(:, 1:2*n), 'fro') / ...
    max(norm(jI, 'fro'), eps);
sensitivityClosure = abs(real(lambda) - (s6 + sE + sI));
slopeProductClosure = abs(slopes.L6Inhibition - ...
    (-slopes.L6Excitation) * slopes.ExcitationInhibition);

storedSlope = NaN;
storedResult = fullfile(outputRoot, 'L6_and_Inhibition_offset', ...
    'full_grid_results.mat');
if isfile(storedResult)
    stored = load(storedResult, 'theorySlope');
    if isfield(stored, 'theorySlope'); storedSlope = stored.theorySlope; end
end

audit = table(real(lambda), imag(lambda), s6, sE, sI, fd6, fdE, fdI, ...
    abs(fd6-s6)/max(abs(s6),eps), ...
    abs(fdE-sE)/max(abs(sE),eps), ...
    abs(fdI-sI)/max(abs(sI),eps), ...
    slopes.L6Inhibition, slopes.L6Excitation, ...
    slopes.ExcitationInhibition, storedSlope, ...
    abs(storedSlope-slopes.L6Inhibition), reconstructionError, ...
    l6ISourceFraction, excitationISourceFraction, ...
    inhibitionNonISourceFraction, sensitivityClosure, ...
    slopeProductClosure, step, ...
    'VariableNames', {'baselineLeadingReal','baselineLeadingImag', ...
    'sensitivityL6','sensitivityExcitation','sensitivityInhibition', ...
    'finiteDifferenceL6','finiteDifferenceExcitation', ...
    'finiteDifferenceInhibition','relativeErrorL6', ...
    'relativeErrorExcitation','relativeErrorInhibition', ...
    'slopeL6Inhibition','slopeL6Excitation', ...
    'slopeExcitationInhibition','storedSlopeL6Inhibition', ...
    'storedSlopeAbsoluteDifference','jacobianReconstructionRelativeError', ...
    'L6ISourceFraction','excitationISourceFraction', ...
    'inhibitionNonISourceFraction','sensitivityClosureAbsoluteError', ...
    'slopeProductClosureAbsoluteError','finiteDifferenceStep'});

writetable(audit, fullfile(outputRoot, 'theoretical_fpp_slope_audit.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(outputRoot, 'theoretical_fpp_slope_audit.mat'), ...
    'audit', 'slopes', '-v7.3');
fprintf(['FPP slopes: L6-I %.12g, L6-E %.12g, E-I %.12g; ' ...
    'stored L6-I difference %.3e.\n'], slopes.L6Inhibition, ...
    slopes.L6Excitation, slopes.ExcitationInhibition, ...
    abs(storedSlope-slopes.L6Inhibition));
end

function sensitivity = local_central_sensitivity(jBaseline, component, step, cfg)
alphaPlus = real(local_leading_mode(jBaseline + step*component, cfg));
alphaMinus = real(local_leading_mode(jBaseline - step*component, cfg));
sensitivity = (alphaPlus-alphaMinus)/(2*step);
end

function [lambda, left, right] = local_leading_mode(matrix, cfg)
count = 6;
options = struct('tol', cfg.EigsTolerance, 'maxit', cfg.EigsMaxIterations, ...
    'p', min(max(cfg.EigsSubspaceDimension, 2*count+8), size(matrix,1)), ...
    'disp', 0, 'isreal', true);
try
    [rightPool, rightValues] = eigs(matrix, count, 'largestreal', options);
catch
    [rightPool, rightValues] = eigs(matrix, count, 'lr', options);
end
rightLambda = diag(rightValues);
[~, index] = max(real(rightLambda));
lambda = rightLambda(index);
right = rightPool(:,index);
right = right / norm(right);
try
    [leftPool, leftValues] = eigs(matrix', count, 'largestreal', options);
catch
    [leftPool, leftValues] = eigs(matrix', count, 'lr', options);
end
leftLambda = diag(leftValues);
[~, leftIndex] = min(abs(leftLambda-conj(lambda)));
left = leftPool(:,leftIndex);
overlap = left' * right;
if abs(overlap) < 1e-12
    error('FppSlope:Eigenvectors', 'Leading left/right overlap is singular.');
end
left = left / conj(overlap);
end
