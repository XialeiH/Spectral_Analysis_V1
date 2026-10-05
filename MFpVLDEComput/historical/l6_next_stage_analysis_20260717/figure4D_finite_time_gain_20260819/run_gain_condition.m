function run_gain_condition(taskId, geometryRoot, outputRoot, codeRoot)
% Compute one Section 3.5 finite-time gain curve for one FPP pathway freeze.

contrasts = [19 42 66 100];
settingNames = {'full','l6_frozen','inhibition_frozen','recurrent_E_frozen'};
[contrastIndex, settingIndex] = ind2sub([numel(contrasts), numel(settingNames)], taskId);

addpath(codeRoot);
cfg = l6ns_config();
cfg.GeometryRoot = geometryRoot;
cfg.Angle = 0;
cfg.Contrast = contrasts(contrastIndex);
cfg.TransientPowerIterations = 24;
cfg.TransientRefinePeak = true;

data = l6ns_load_endpoints(cfg);
pathway = l6ns_two_pathway_data(data);
switch settingNames{settingIndex}
    case 'full'
        J = pathway.JBaseline;
    case 'l6_frozen'
        J = pathway.JRest + pathway.JI;
    case 'inhibition_frozen'
        J = pathway.JRest + pathway.J6;
    case 'recurrent_E_frozen'
        J = pathway.J6 + pathway.JI;
    otherwise
        error('Unknown pathway setting.');
end

n = data.PopulationSize;
fixedPoint = data.FixedPoint(:);
populationRms = [norm(fixedPoint(1:n)), ...
    norm(fixedPoint(n+(1:n))), norm(fixedPoint(2*n+(1:n)))] / sqrt(n);
populationRms = max(populationRms, 1e-8);
weightVector = [ones(n,1)/populationRms(1); ...
    ones(n,1)/populationRms(2); ones(n,1)/populationRms(3)];

timeMs = unique([0:2.5:20 25:5:60 70:10:120]);
gainResult = l6ns_continuous_transient(J, timeMs, cfg, weightVector, cfg.TauMs);

metadata = struct();
metadata.Section35Definition = 'G(t)=Wout^(1/2)*Cobs*exp(L*t)*Bpert*Win^(-1/2)';
metadata.Specialization = ['Bpert=Cobs=I; Wout^(1/2)=Win^(1/2)=' ...
    'diag(1./populationRms repeated by population)'];
metadata.Generator = 'L=(J-I)/tau';
metadata.TauMs = cfg.TauMs;
metadata.AngleDeg = cfg.Angle;
metadata.Contrast = cfg.Contrast;
metadata.Setting = settingNames{settingIndex};
metadata.PopulationRms = populationRms;
metadata.JacobianDecomposition = ...
    'Jfull=Jrest+J6+JI; freezes remove J6, JI, or Jrest at the same fixed point';
metadata.PathwayReconstructionError = pathway.ReconstructionError;

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
outputFile = fullfile(outputRoot, sprintf('gain_contr%03d_%s.mat', ...
    cfg.Contrast, settingNames{settingIndex}));
save(outputFile, 'gainResult', 'metadata', '-v7.3');
fprintf('Saved %s\n', outputFile);
end
