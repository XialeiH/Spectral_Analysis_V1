function paths = mechanism_initialize()
% Resolve project paths and reproduce the Torch helper precedence locally.

paths.MechanismRoot = fileparts(mfilename('fullpath'));
paths.RealTuningRoot = fileparts(paths.MechanismRoot);
paths.L6AnalysisRoot = fileparts(paths.RealTuningRoot);
paths.InsertRoot = fileparts(paths.L6AnalysisRoot);
paths.MainRoot = fullfile(paths.InsertRoot,'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
paths.RuntimeRoot = fullfile(paths.MechanismRoot,'runtime_h96');
paths.SetupFile = fullfile(paths.L6AnalysisRoot, ...
    'results_global_bifurcation_20260722','global_bifurcation_setup.mat');
paths.ResultsRoot = fullfile(paths.RealTuningRoot,'results_full_20260722_055032');
paths.OutputRoot = fullfile(paths.MechanismRoot,'results');
paths.FigureRoot = fullfile(paths.OutputRoot,'figures');
if ~exist(paths.OutputRoot,'dir'); mkdir(paths.OutputRoot); end
if ~exist(paths.FigureRoot,'dir'); mkdir(paths.FigureRoot); end

addpath(paths.MainRoot);
addpath(fullfile(paths.MainRoot,'Utils'));
addpath(paths.RuntimeRoot,'-begin');
addpath(paths.L6AnalysisRoot,'-begin');
addpath(paths.RealTuningRoot,'-begin');
addpath(paths.MechanismRoot,'-begin');
end
