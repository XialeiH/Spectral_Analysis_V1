function prepare_Figure_5E_2_tuning
% Recompute Figure 5E.2 tuning curves with the current project pixel remap.

artifactRoot = fileparts(mfilename('fullpath'));
projectRoot = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis/matlab-inserting_into_CG_model'];
mainRoot = fullfile(projectRoot, ...
    'Complete_Code_for_Paper3', 'NYU-Vision-2Drive-main');
analysisRoot = fullfile(projectRoot, 'l6_next_stage_analysis_20260717');
figure60Root = fullfile(analysisRoot, ...
    'figure6_0_multimetric_20260807');
fixedPointRoot = fullfile(analysisRoot, ...
    'top_singular_hcnorm_grids_20260803');
h96Runtime = fullfile(analysisRoot, 'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722', 'runtime_h96');
dataFile = fullfile(mainRoot, 'Figures', ...
    'spectral_analysis_eigenvalue_eigenvectors', 'L6 and Inhibition', ...
    'L6_and_Inhibition_offset', 'Figure6_0_FourContrast_FiveMetric', ...
    'figure6_0_contrast100_dataset.mat');
outputFile = fullfile(artifactRoot, 'Figure_5E_2_corrected_tuning.mat');

addpath(mainRoot);
addpath(fullfile(mainRoot, 'Utils'));
addpath(h96Runtime, '-begin');
addpath(analysisRoot, '-begin');
addpath(fixedPointRoot, '-begin');
addpath(figure60Root, '-begin');

loaded = load(dataFile, 'setup');
setup = loaded.setup;
actualBeta6 = [0, 0.072, 0.296];
actualBetaI = [0, 0.040, 0.168];
canonicalStates = zeros([size(setup.BaselineCanonicalStates), 3]);
canonicalStates(:, :, 1) = double(setup.BaselineCanonicalStates);
canonicalResiduals = zeros(numel(setup.CanonicalAngles), 3);
canonicalIterations = zeros(numel(setup.CanonicalAngles), 3);

for conditionIndex = 2:3
    context = setup.Context;
    for angleIndex = 1:numel(setup.CanonicalAngles)
        angleContext = context;
        angleContext.OrientationUse = setup.CanonicalAngles(angleIndex) * ...
            ones(size(context.OrientationUse));
        phi = @(state) l6ns_phi(state, -actualBeta6(conditionIndex), ...
            angleContext, [1 1], 1 + actualBetaI(conditionIndex), 'true');
        fixed = real_tuning_fixed_point_tolerance(phi, ...
            setup.BaselineCanonicalStates(:, angleIndex), ...
            angleContext.RelaxationP, 1e-4);
        if ~fixed.Converged
            error('Figure5E2:TuningFixedPoint', ...
                ['Condition %d, canonical angle %.2f deg failed with ' ...
                'residual %.6g (%s).'], conditionIndex, ...
                setup.CanonicalAngles(angleIndex), fixed.Residual, ...
                fixed.Termination);
        end
        canonicalStates(:, angleIndex, conditionIndex) = fixed.State;
        canonicalResiduals(angleIndex, conditionIndex) = fixed.Residual;
        canonicalIterations(angleIndex, conditionIndex) = fixed.Iterations;
    end
end

tuningCurves = zeros(numel(setup.FullAngles), 2, 3);
for conditionIndex = 1:3
    tuningCurves(:, :, conditionIndex) = figure6_0_tuning_curves( ...
        canonicalStates(:, :, conditionIndex), ...
        setup.CanonicalAngles, setup.CWeight);
end

pixelRows = [5, 1];
pixelCols = [10, 10];
preferredAnglesDeg = [0, 22.5];
fullAngles = double(setup.FullAngles(:));
save(outputFile, 'tuningCurves', 'fullAngles', 'pixelRows', ...
    'pixelCols', 'preferredAnglesDeg', 'actualBeta6', 'actualBetaI', ...
    'canonicalResiduals', 'canonicalIterations', '-v7.3');
fprintf('Saved corrected Figure 5E.2 tuning curves to %s.\n', outputFile);
fprintf('Maximum fixed-point residual: %.6g.\n', ...
    max(canonicalResiduals, [], 'all'));
end
