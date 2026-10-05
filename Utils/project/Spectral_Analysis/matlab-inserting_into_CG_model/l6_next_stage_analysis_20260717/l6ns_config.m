function cfg = l6ns_config()
% Configuration shared by the five L6 next-stage analyses.

thisDir = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');

cfg = struct();
cfg.CodeRoot = thisDir;
cfg.GeometryRoot = getenv('L6NS_GEOMETRY_ROOT');
if isempty(cfg.GeometryRoot)
    cfg.GeometryRoot = fullfile(thisDir, 'geometry_mat');
end
cfg.OutputRoot = getenv('L6NS_OUTPUT_ROOT');
if isempty(cfg.OutputRoot)
    cfg.OutputRoot = fullfile(thisDir, 'results');
end

cfg.Angle = local_env_number('L6NS_ANGLE', 0);
cfg.Contrast = local_env_number('L6NS_CONTRAST', 100);
cfg.CriticalBoundary = 1;
cfg.StabilityConvention = 'continuous_phi_minus_identity';
cfg.StateLabels = {'S','C','I'};
cfg.ExcitatoryCWeight = 0.3077;
cfg.LifPeakTimeMs = 23;
cfg.LifPsthBinMs = 2;
cfg.LifPeakWindow = 12;
cfg.LifReferenceFile = 'DriveWkSp_Ang0.0_total200ms1s.mat';
cfg.TauMs = 10.3402405296839;
cfg.PhysicalTimeGridMs = unique([0:2.5:20, 25:5:60, 70:10:200]);
cfg.GridTimeGridMs = unique([0:5:30, 40:10:120]);
cfg.Gamma6Grid = 0:0.1:1.5;
cfg.GammaIGrid = 0.5:0.05:1.5;
cfg.InteractionGamma6Offsets = [-0.2 -0.1 0 0.1 0.2];
cfg.MatchedStabilityCount = 3;
cfg.ContrastProbeStep = 1;
cfg.OrientationProbeStepDeg = 0.25;
cfg.SurroundProbeStep = 1;
cfg.TransientPowerIterations = 40;
cfg.TransientRefinePeak = true;
cfg.TransientSvdsTolerance = 1e-7;
cfg.TransientSvdsMaxIterations = 100;
cfg.TransientSvdsSubspaceDimension = 16;

cfg.DenseWGrid = -0.35:0.005:0.05;
cfg.AsymptoticWGrid = [-1 -2 -5 -10 -20 -50];
cfg.SubspaceWGrid = [-0.30 -0.20 -0.16 -0.10 0 0.5 0.8 1 2 10 20];
cfg.SubspaceDimensions = [4 8 20 50];
cfg.BranchPoolSize = 14;
cfg.PositiveSubspaceSize = 20;
cfg.ActiveSingularSize = 40;
cfg.ComponentAtlasSize = 24;
cfg.AnalysisModeCount = 100;
cfg.SingularModeCount = 100;
cfg.BaseModalCount = 80;
cfg.ReducedDimensions = [8 16 32 64 96];
cfg.TransientTimeGrid = linspace(0, 40, 41);
cfg.ComputeFullSingularSpectrum = true;
cfg.ComputeFullEigenSpectrum = true;

cfg.EigsTolerance = 1e-9;
cfg.EigsMaxIterations = 1600;
cfg.EigsSubspaceDimension = 80;
cfg.RandomSeed = 7;

cfg.Smoke = local_env_boolean('L6NS_SMOKE', false);
if cfg.Smoke
    cfg.DenseWGrid = [-0.30 -0.22 -0.18 -0.16 -0.14 -0.10 0 0.05];
    cfg.AsymptoticWGrid = [-1 -5 -20];
    cfg.SubspaceWGrid = [-0.20 -0.16 0 1];
    cfg.SubspaceDimensions = [4 8];
    cfg.BranchPoolSize = 8;
    cfg.PositiveSubspaceSize = 8;
    cfg.ActiveSingularSize = 10;
    cfg.ComponentAtlasSize = 6;
    cfg.AnalysisModeCount = 16;
    cfg.SingularModeCount = 16;
    cfg.BaseModalCount = 16;
    cfg.ReducedDimensions = [8 16];
    cfg.TransientTimeGrid = [0 1 2 4 8];
    cfg.PhysicalTimeGridMs = [0 5 10 20 40];
    cfg.GridTimeGridMs = [0 10 20 40];
    cfg.Gamma6Grid = [0 0.5 1 1.25 1.5];
    cfg.GammaIGrid = [0.5 0.8 1 1.2 1.5];
    cfg.InteractionGamma6Offsets = [-0.1 -0.05 0 0.05 0.1];
    cfg.MatchedStabilityCount = 2;
    cfg.TransientPowerIterations = 20;
    cfg.ComputeFullSingularSpectrum = false;
    cfg.ComputeFullEigenSpectrum = false;
end

if ~exist(cfg.OutputRoot, 'dir')
    mkdir(cfg.OutputRoot);
end

function value = local_env_boolean(name, defaultValue)
raw = strtrim(lower(getenv(name)));
if isempty(raw)
    value = defaultValue;
elseif any(strcmp(raw, {'1','true','yes','on'}))
    value = true;
elseif any(strcmp(raw, {'0','false','no','off'}))
    value = false;
else
    error('%s must be a boolean value.', name);
end
end
end

function value = local_env_number(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
    return
end
value = str2double(raw);
if ~isfinite(value)
    error('%s must contain a finite number.', name);
end
end
