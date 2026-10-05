function summary = run_jrest_top4_atlas(sourceFile, outputRoot)
% Plot the top-four positive-edge modes of J_rest at contrast 100, angle 0.

if nargin < 1 || isempty(sourceFile)
    sourceFile = getenv('JREST_SOURCE_MAT');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = getenv('JREST_OUTPUT_ROOT');
end
if isempty(sourceFile) || ~isfile(sourceFile)
    error('JRestAtlas:Source', ...
        'JREST_SOURCE_MAT must identify the audited two_pathway_result.mat.');
end
if isempty(outputRoot)
    error('JRestAtlas:Output', 'JREST_OUTPUT_ROOT is required.');
end

loaded = load(sourceFile, 'result');
if ~isfield(loaded, 'result') || ~isfield(loaded.result, 'Pathway')
    error('JRestAtlas:Pathway', 'The source file does not contain result.Pathway.');
end
pathway = loaded.result.Pathway;
required = {'JRest','J6','JI','JBaseline'};
for index = 1:numel(required)
    if ~isfield(pathway, required{index})
        error('JRestAtlas:MissingField', ...
            'result.Pathway.%s is missing.', required{index});
    end
end

jRest = sparse(pathway.JRest);
reconstructed = jRest + sparse(pathway.J6) + sparse(pathway.JI);
reconstructionRelativeError = norm(reconstructed - sparse(pathway.JBaseline), 'fro') / ...
    max(norm(sparse(pathway.JBaseline), 'fro'), eps);
if reconstructionRelativeError >= 1e-13
    error('JRestAtlas:Reconstruction', ...
        'JRest + J6 + JI does not reconstruct JBaseline: relative error %.3e.', ...
        reconstructionRelativeError);
end
clear reconstructed pathway loaded

cfg = l6ns_config();
cfg.EigsTolerance = 1e-10;
cfg.EigsMaxIterations = 2200;
cfg.EigsSubspaceDimension = 80;
rng(3100, 'twister');
modes = l6ns_eigenpairs(jRest, 4, 'largestreal', cfg);

angleValue = 0;
contrast = 100;
n = size(modes.Right, 1) / 3;
side = round(sqrt(n));
if n ~= round(n) || side * side ~= n
    error('JRestAtlas:MapSize', ...
        'The three-population state does not form square maps.');
end
mapSize = [side side];
populationLabels = {'S','C','I','E = 0.6923 S + 0.3077 C'};

figureRoot = fullfile(outputRoot, 'figures');
dataRoot = fullfile(outputRoot, 'data');
if ~exist(figureRoot, 'dir'); mkdir(figureRoot); end
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [40 40 2200 1850]);
layout = tiledlayout(figureHandle, 4, 4, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
colormap(figureHandle, parula(256));

energyFractions = zeros(4, 3);
for modeIndex = 1:4
    vector = modes.Right(:, modeIndex);
    s = vector(1:n);
    c = vector(n + (1:n));
    i = vector(2*n + (1:n));
    e = 0.6923 * s + 0.3077 * c;
    maps = {reshape(real(s), mapSize), reshape(real(c), mapSize), ...
        reshape(real(i), mapSize), reshape(real(e), mapSize)};
    limit = max(cellfun(@(map) max(abs(map(:))), maps));
    if limit <= eps; limit = 1; end

    energy = [sum(abs(s).^2), sum(abs(c).^2), sum(abs(i).^2)];
    energyFractions(modeIndex, :) = energy / max(sum(energy), eps);
    for populationIndex = 1:4
        axisHandle = nexttile(layout);
        imagesc(axisHandle, maps{populationIndex});
        axis(axisHandle, 'image');
        axis(axisHandle, 'off');
        clim(axisHandle, [-limit limit]);
        colorbar(axisHandle, 'FontSize', 9);
        if modeIndex == 1
            title(axisHandle, populationLabels{populationIndex}, ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
        if populationIndex == 1
            text(axisHandle, -0.08, 0.5, ...
                sprintf('mode %d: %.5g%+.5gi', modeIndex, ...
                real(modes.Lambda(modeIndex)), imag(modes.Lambda(modeIndex))), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', 'Rotation', 90, ...
                'Clipping', 'off', 'FontSize', 10);
        end
    end
end

title(layout, { ...
    'J_rest positive-edge top 4 | angle 0.00 deg | contrast 100', ...
    sprintf(['J_rest = J_baseline - J_6 - J_I | dynamic L6 and inhibitory ' ...
    'feedback removed; S/C/I states retained | max Re(lambda)=%.9f'], ...
    real(modes.Lambda(1)))}, ...
    'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');
stem = 'Top4_J_rest_Angle_0p00deg_Contrast_100';
l6ns_save_figure(figureHandle, figureRoot, stem);
close(figureHandle);

mode = (1:4).';
lambdaReal = real(modes.Lambda(:));
lambdaImag = imag(modes.Lambda(:));
conditionNumber = modes.ConditionNumber(:);
rightResidual = modes.RightResidual(:);
leftResidual = modes.LeftResidual(:);
S_energyFraction = energyFractions(:, 1);
C_energyFraction = energyFractions(:, 2);
I_energyFraction = energyFractions(:, 3);
summary = table(mode, lambdaReal, lambdaImag, conditionNumber, ...
    rightResidual, leftResidual, S_energyFraction, C_energyFraction, ...
    I_energyFraction);
writetable(summary, fullfile(dataRoot, ...
    'top4_J_rest_summary_angle0p00_contrast100.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

metadata = struct( ...
    'Definition', 'JRest=JBaseline-J6-JI', ...
    'Interpretation', ...
        ['Dynamic L6 and inhibitory derivative-feedback blocks are removed; ' ...
         'the S, C, and I state populations remain in the Jacobian.'], ...
    'Angle', angleValue, ...
    'Contrast', contrast, ...
    'ReconstructionRelativeError', reconstructionRelativeError, ...
    'SavedFullMatrix', false, ...
    'SourceFile', sourceFile);
modeData = struct('Lambda', modes.Lambda, 'Right', modes.Right, ...
    'Left', modes.Left, 'ConditionNumber', modes.ConditionNumber, ...
    'RightResidual', modes.RightResidual, 'LeftResidual', modes.LeftResidual);
save(fullfile(dataRoot, 'top4_J_rest_data_angle0p00_contrast100.mat'), ...
    'summary', 'metadata', 'modeData', '-v7.3');

fprintf(['Saved J_rest top-four atlas. max Re(lambda)=%.9f, ' ...
    'reconstruction relative error %.3e.\n'], ...
    real(modes.Lambda(1)), reconstructionRelativeError);
end
