function prepare_Figure_3G_data()
% Compute the full 4,800-state alpha-by-eta intervention grid.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
modelRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model');
analysisRoot = fullfile(modelRoot, 'l6_next_stage_analysis_20260717');
controlRoot = fullfile(analysisRoot, ...
    'spatial_symmetry_zero_controls_20260811');
responseRuntime = fullfile(analysisRoot, 'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722', 'runtime_h96');
gradientRuntime = fullfile(modelRoot, 'geometry_analysis', ...
    'torch_h96_dedupe_array', 'Utils');
utilsRoot = fullfile(modelRoot, 'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main', 'Utils');
outputRoot = fullfile(modelRoot, 'report_figure_artifacts', ...
    '3_Weighted_Jacobian', 'Figure_3G');
cacheRoot = fullfile(outputRoot, 'spectrum_cache');
if ~isfolder(cacheRoot); mkdir(cacheRoot); end

addpath(utilsRoot);
addpath(analysisRoot);
addpath(controlRoot);
addpath(gradientRuntime);
addpath(responseRuntime, '-begin');

loaded = load(fullfile(analysisRoot, ...
    'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat'), 'setup');
setup = loaded.setup;
alphaValues = 0:0.1:0.5;
etaValues = 0:0.2:1;
threshold = 0.05;
nearZeroCounts = nan(numel(etaValues), numel(alphaValues));
maxRealEigenvalue = nan(size(nearZeroCounts));

for alphaIndex = 1:numel(alphaValues)
    alpha = alphaValues(alphaIndex);
    fprintf('Constructing fixed-gain Gaussian control alpha = %.1f.\n', alpha);
    Jalpha = sparse(l6ns_spatial_smoothing_control_jacobian( ...
        setup, alpha, alpha));
    if alphaIndex == 1
        relativeError = norm(Jalpha - setup.Pathway.JBaseline, 'fro') / ...
            norm(setup.Pathway.JBaseline, 'fro');
        fprintf('alpha=0 relative Jacobian error: %.3e.\n', relativeError);
        assert(relativeError < 1e-9, ...
            'alpha=0 does not reproduce the baseline Jacobian.');
    end

    n = size(Jalpha, 1) / 3;
    e = 1:(2*n);
    i = (2*n+1):(3*n);
    for etaIndex = 1:numel(etaValues)
        eta = etaValues(etaIndex);
        cacheFile = fullfile(cacheRoot, sprintf( ...
            'spectrum_alpha_%s_eta_%s.mat', tag(alpha), tag(eta)));
        if isfile(cacheFile)
            cached = load(cacheFile, 'eigenvalues');
            eigenvalues = cached.eigenvalues;
            fprintf('Loaded alpha %.1f, eta %.1f.\n', alpha, eta);
        else
            J = Jalpha;
            crossScale = sqrt(eta);
            J(e,i) = crossScale * Jalpha(e,i);
            J(i,e) = crossScale * Jalpha(i,e);
            fprintf('Computing alpha %.1f, eta %.1f full spectrum.\n', ...
                alpha, eta);
            timer = tic;
            eigenvalues = eig(full(J), 'vector');
            elapsedSeconds = toc(timer);
            save(cacheFile, 'alpha', 'eta', 'eigenvalues', ...
                'elapsedSeconds', '-v7.3');
            fprintf('Completed in %.1f s.\n', elapsedSeconds);
        end
        assert(numel(eigenvalues) == 4800, ...
            'Unexpected eigenspectrum dimension.');
        nearZeroCounts(etaIndex, alphaIndex) = ...
            nnz(abs(eigenvalues) <= threshold);
        maxRealEigenvalue(etaIndex, alphaIndex) = max(real(eigenvalues));
        fprintf('N_0.05 = %d; max Re(lambda) = %.6f.\n', ...
            nearZeroCounts(etaIndex, alphaIndex), ...
            maxRealEigenvalue(etaIndex, alphaIndex));
        clear eigenvalues J
    end
    clear Jalpha
end

baselineCount = nearZeroCounts(etaValues == 1, alphaValues == 0);
noLoopCount = nearZeroCounts(etaValues == 0, alphaValues == 0);
flatCount = nearZeroCounts(etaValues == 1, alphaValues == 0.5);
assert(baselineCount == 2976, ...
    'Baseline near-zero count changed: %d.', baselineCount);
deltaEI = baselineCount - noLoopCount;
deltaGaussian = baselineCount - flatCount;
outputFile = fullfile(outputRoot, 'figure_3g_data.mat');
save(outputFile, 'alphaValues', 'etaValues', 'threshold', ...
    'nearZeroCounts', 'maxRealEigenvalue', 'baselineCount', ...
    'noLoopCount', 'flatCount', 'deltaEI', 'deltaGaussian', ...
    'cacheRoot', '-v7.3');
fprintf('Saved %s\n', outputFile);
fprintf('Delta_EI = %d (%.4f baseline); Delta_Gaussian = %d (%.4f baseline).\n', ...
    deltaEI, deltaEI/baselineCount, ...
    deltaGaussian, deltaGaussian/baselineCount);
end

function value = tag(number)
value = strrep(sprintf('%.1f', number), '.', 'p');
end
