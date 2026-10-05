function run_spectral_edge_condition(taskId, geometryRoot, outputRoot, codeRoot)
% Compute the four FPP pathway-freeze spectral edges for one condition.

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[angleIndex, contrastIndex] = ind2sub([numel(angles), numel(contrasts)], taskId);

addpath(codeRoot);
cfg = l6ns_config();
cfg.GeometryRoot = geometryRoot;
cfg.Angle = angles(angleIndex);
cfg.Contrast = contrasts(contrastIndex);

data = l6ns_load_endpoints(cfg);
pathway = l6ns_two_pathway_data(data);
settingNames = {'full','l6_frozen','inhibition_frozen','recurrent_E_frozen'};
jacobians = {pathway.JBaseline, pathway.JRest + pathway.JI, ...
    pathway.JRest + pathway.J6, pathway.J6 + pathway.JI};

maxRealLambda = nan(1, numel(jacobians));
eigenResidual = nan(1, numel(jacobians));
for settingIndex = 1:numel(jacobians)
    [maxRealLambda(settingIndex), eigenResidual(settingIndex)] = ...
        local_spectral_edge(jacobians{settingIndex}, cfg);
end

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
result = table(repmat(taskId, numel(settingNames), 1), ...
    repmat(cfg.Angle, numel(settingNames), 1), ...
    repmat(cfg.Contrast, numel(settingNames), 1), string(settingNames(:)), ...
    maxRealLambda(:), eigenResidual(:), ...
    repmat(pathway.ReconstructionError, numel(settingNames), 1), ...
    'VariableNames', {'TaskId','AngleDeg','Contrast','Setting', ...
    'MaxRealLambda','EigenResidual','PathwayReconstructionError'});
outputFile = fullfile(outputRoot, sprintf('spectral_edge_task_%02d.tsv', taskId));
writetable(result, outputFile, 'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved %s\n', outputFile);
end

function [value, residual] = local_spectral_edge(jacobian, cfg)
opts = struct('tol', cfg.EigsTolerance, 'maxit', cfg.EigsMaxIterations, ...
    'p', min(max(cfg.EigsSubspaceDimension, 32), size(jacobian, 1)), ...
    'disp', 0);
[vectors, values, flag] = eigs(jacobian, 6, 'largestreal', opts);
if flag ~= 0
    error('eigs did not converge (flag %d).', flag);
end
lambda = diag(values);
[value, index] = max(real(lambda));
vector = vectors(:, index);
residual = norm(jacobian * vector - lambda(index) * vector) / ...
    max((norm(jacobian, 'fro') + abs(lambda(index))) * norm(vector), eps);
end
