function evaluate_h96_paper3_snn_pairs()
% Evaluate Paper3 SNN pixel inputs with the corrected h96 NN surrogate.

runtimeDir = getenv('H96_CORRECTED_RUNTIME');
if isempty(runtimeDir)
    runtimeDir = ['/scratch/xh2906/librarySCI_runs/' ...
        'sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/' ...
        'hpc_runtime_bundle_extended_corrected_h96_ikp7_l6p5_matchedlib_20260702_010052'];
end
utilsDir = getenv('H96_CURRENT_UTILS');
simulationDir = getenv('SNN_OUTPUT_DIR');
outputFile = getenv('PREMODEL_PAIR_DATA');
if isempty(utilsDir) || isempty(simulationDir) || isempty(outputFile)
    error('H96_CURRENT_UTILS, SNN_OUTPUT_DIR, and PREMODEL_PAIR_DATA are required.');
end
addpath(utilsDir, '-begin');
addpath(runtimeDir, '-begin');

expectedHash = '41c7fdc837617cd66ee81b0a7f7dc2bc3af45ed8c38e67bc3e85abed8b6bb161';
resolvedResponse = which('LocalResponse_6D_MLP_prefAngle');
assert(startsWith(resolvedResponse, runtimeDir), ...
    'Wrong h96 LocalResponse resolved: %s', resolvedResponse);
assertSha256(resolvedResponse, expectedHash);
assertSha256(which('predict_pref6D_S'), ...
    '0eaf440aebc1ad3194c01efb3a807ab8eb03d8dff223a64fdd1565935da6049a');
assertSha256(which('predict_pref6D_C'), ...
    '1efcf8260cae702885f775e8875e098a9658d09aad7f74631ff701ec28709ae2');
assertSha256(which('predict_pref6D_I'), ...
    '6be47fb16e2b913f6cb7957b6b9daad9b6862f9b85fdd05822c4b54c9ac7c971');
fprintf('Corrected h96 runtime: %s\n', resolvedResponse);
wC = 0.3077;
angles = [0, 7.5, 15, 22.5];
baseReplicate = 1;
extraReplicate = 2;
extraCountPerAngle = 270;
evenIndices = unique(round(linspace(1, 900, extraCountPerAngle)));
assert(numel(evenIndices) == extraCountPerAngle, 'Even 30%% sample has wrong size.');

E_premodel = []; E_h96 = []; I_premodel = []; I_h96 = [];
angle = []; pixelIndex = []; replicate = []; isExtra30 = [];
inputRanges = zeros(numel(angles) * 2, 18);
rangeRow = 0;
for ai = 1:numel(angles)
    for rep = [baseReplicate, extraReplicate]
        file = fullfile(simulationDir, sprintf( ...
            'Paper3NWSimulationPix_%.1fdeg_rep%02d_h96_validation.mat', angles(ai), rep));
        assert(isfile(file), 'Missing SNN simulation: %s', file);
        D = load(file, 'NWSmlt', 'metadata');
        [networkE, predictedE, networkI, predictedI] = ...
            evaluateOne(D.NWSmlt, angles(ai), wC);
        rangeRow = rangeRow + 1;
        inputRanges(rangeRow,:) = populationInputRanges(D.NWSmlt);
        if rep == baseReplicate
            keep = (1:900)';
            extraFlag = false(900,1);
        else
            keep = evenIndices(:);
            extraFlag = true(numel(keep),1);
        end
        E_premodel = [E_premodel; networkE(keep)]; %#ok<AGROW>
        E_h96 = [E_h96; predictedE(keep)]; %#ok<AGROW>
        I_premodel = [I_premodel; networkI(keep)]; %#ok<AGROW>
        I_h96 = [I_h96; predictedI(keep)]; %#ok<AGROW>
        angle = [angle; repmat(angles(ai),numel(keep),1)]; %#ok<AGROW>
        pixelIndex = [pixelIndex; keep]; %#ok<AGROW>
        replicate = [replicate; repmat(rep,numel(keep),1)]; %#ok<AGROW>
        isExtra30 = [isExtra30; extraFlag]; %#ok<AGROW>
    end
end

assert(numel(E_premodel) == 4680, 'Expected 4680 points after adding 30%%.');
assert(all(isfinite([E_premodel;E_h96;I_premodel;I_h96])), 'Nonfinite paired data.');
paired = struct('E_premodel',E_premodel,'E_h96',E_h96, ...
    'I_premodel',I_premodel,'I_h96',I_h96,'angleDeg',angle, ...
    'pixelIndex',pixelIndex,'replicate',replicate,'isExtra30',isExtra30);
provenance = struct('runtimeDir',runtimeDir,'localResponse',resolvedResponse, ...
    'localResponseSha256',expectedHash,'utilsDir',utilsDir, ...
    'simulationDir',simulationDir,'anglesDeg',angles, ...
    'baseReplicate',baseReplicate,'extraReplicate',extraReplicate, ...
    'basePointCount',3600,'extraPointCount',1080,'extraFraction',0.30, ...
    'extraPixelIndices',evenIndices,'inputRanges',inputRanges, ...
    'comparison','raw Paper3 SNN firing rate versus raw h96 local-response output', ...
    'postPredictionSaturationApplied',false);
save(outputFile, 'paired', 'provenance', '-v7.3');
fprintf('Saved %s with %d paired points per population.\n', outputFile, numel(E_premodel));
end

function [networkE,predictedE,networkI,predictedI] = evaluateOne(NWSmlt, angleDeg, wC)
networkE = NWSmlt.FS * (1-wC) + NWSmlt.FC * wC;
networkI = NWSmlt.FI;
alpha = angleDeg * ones(900,1);
contrast = 100 * ones(900,1);
rawS = evalPopulation('S',alpha,[NWSmlt.LGNWeightsS,zeros(900,1)], ...
    NWSmlt.L4EByPopulation(:,1),NWSmlt.L4IByPopulation(:,1), ...
    NWSmlt.L6ByPopulation(:,1),contrast);
rawC = evalPopulation('C',alpha,[NWSmlt.LGNWeightsC,zeros(900,1)], ...
    NWSmlt.L4EByPopulation(:,2),NWSmlt.L4IByPopulation(:,2), ...
    NWSmlt.L6ByPopulation(:,2),contrast);
rawI = evalPopulation('I',alpha,[NWSmlt.LGNWeightsI,zeros(900,1)], ...
    NWSmlt.L4EByPopulation(:,3),NWSmlt.L4IByPopulation(:,3), ...
    NWSmlt.L6ByPopulation(:,3),contrast);
predictedE = rawS*(1-wC) + rawC*wC;
predictedI = rawI;
end

function result = evalPopulation(cellType,alpha,pixLGN,l4E,l4I,l6,contrast)
allCategories = zeros(numel(alpha),5);
for category = 1:5
    allCategories(:,category) = LocalResponse_6D_MLP_prefAngle(cellType,alpha, ...
        category*ones(size(alpha)),l4E,l4I,l6,contrast);
end
result = sum(pixLGN .* allCategories, 2);
end

function ranges = populationInputRanges(NWSmlt)
ranges = zeros(1,18);
offset = 0;
for population = 1:3
    ranges(offset+(1:6)) = [min(NWSmlt.L4EByPopulation(:,population)), ...
        max(NWSmlt.L4EByPopulation(:,population)), ...
        min(NWSmlt.L4IByPopulation(:,population)), ...
        max(NWSmlt.L4IByPopulation(:,population)), ...
        min(NWSmlt.L6ByPopulation(:,population)), ...
        max(NWSmlt.L6ByPopulation(:,population))];
    offset = offset + 6;
end
end

function assertSha256(file, expected)
assert(~isempty(file) && isfile(file), 'Missing required MATLAB helper: %s', file);
[status, hashText] = system(sprintf('sha256sum "%s"', file));
if status ~= 0
    [status, hashText] = system(sprintf('shasum -a 256 "%s"', file));
end
assert(status == 0 && startsWith(strtrim(hashText), expected), ...
    'SHA256 mismatch for %s: %s', file, strtrim(hashText));
end
