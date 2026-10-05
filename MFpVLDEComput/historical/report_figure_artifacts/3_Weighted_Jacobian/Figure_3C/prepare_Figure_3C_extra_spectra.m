function prepare_Figure_3C_extra_spectra()
% Cache complete spectra for the three added Gaussian-kernel controls.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
modelRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model');
analysisRoot = fullfile(modelRoot, 'l6_next_stage_analysis_20260717');
controlRoot = fullfile(analysisRoot, ...
    'spatial_symmetry_zero_controls_20260811');
gkeRoot = fullfile(analysisRoot, 'gaussian_kernel_experiments_20260813');
responseRuntime = fullfile(analysisRoot, 'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722', 'runtime_h96');
gradientRuntime = fullfile(modelRoot, 'geometry_analysis', ...
    'torch_h96_dedupe_array', 'Utils');
utilsRoot = fullfile(modelRoot, 'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main', 'Utils');
cacheRoot = fullfile(modelRoot, 'report_figure_artifacts', ...
    '3_Weighted_Jacobian', 'Figure_3C', 'gke_cache');
equilibriumRoot = fullfile(cacheRoot, 'equilibrium_only', 'data');
spectrumRoot = fullfile(cacheRoot, 'complete_spectra');
if ~exist(spectrumRoot, 'dir'); mkdir(spectrumRoot); end

addpath(utilsRoot);
addpath(analysisRoot);
addpath(controlRoot);
addpath(gkeRoot);
addpath(gradientRuntime);
addpath(responseRuntime, '-begin');

setupData = load(fullfile(analysisRoot, ...
    'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat'), 'setup');
setup = setupData.setup;
names = {'shape_inscribed_triangle', 'shape_left_half', ...
    'kernel_noise_cv0p30'};

for index = 1:numel(names)
    outputFile = fullfile(spectrumRoot, [names{index} '.mat']);
    if isfile(outputFile)
        fprintf('Keeping existing %s.\n', outputFile);
        continue
    end
    equilibriumData = load(fullfile(equilibriumRoot, ...
        [names{index} '.mat']), 'fixedPoint');
    [operators, ~] = gke_build_operators(setup.Context, names{index});
    jacobian = sparse(l6ns_controlled_jacobian( ...
        setup, operators, equilibriumData.fixedPoint));
    fprintf('%s: computing complete 4,800-state eigenspectrum.\n', names{index});
    eigenvalues = eig(full(jacobian), 'vector');
    leadingSingularValue = svds(jacobian, 1, 'largest', ...
        struct('tol', 1e-9, 'maxit', 6000, 'p', 40, 'disp', 0));
    fixedPoint = equilibriumData.fixedPoint;
    save(outputFile, 'fixedPoint', 'eigenvalues', ...
        'leadingSingularValue', '-v7.3');
    fprintf('%s saved: max Re(lambda) %.8f.\n', ...
        names{index}, max(real(eigenvalues)));
    clear jacobian eigenvalues fixedPoint
end
end
