function insetData = build_Figure_5D_inset_data
% Track leading spectral branches and assemble rank-4 angle diagnostics.

thisRoot = fileparts(mfilename('fullpath'));
pointsRoot = fullfile(thisRoot, 'subspace_s0025', 'points');
files = dir(fullfile(pointsRoot, 'subspace_modes_*.mat'));
[~, order] = sort({files.name});
files = files(order);
if numel(files) ~= 22
    error('Figure5D:ModeFiles', 'Expected 22 subspace mode files.');
end

branchCount = 8;
candidateCount = 16;
pointCount = numel(files);
arcLength = zeros(pointCount, 1);
principalAnglesDeg = zeros(pointCount, 4);
lambdaByPoint = cell(pointCount, 1);
vectorsByPoint = cell(pointCount, 1);
for pointIndex = 1:pointCount
    loaded = load(fullfile(files(pointIndex).folder, files(pointIndex).name), ...
        'arcLength', 'lambda', 'vectors', 'principalAnglesRank4Deg');
    arcLength(pointIndex) = loaded.arcLength;
    lambdaByPoint{pointIndex} = loaded.lambda(:);
    vectorsByPoint{pointIndex} = normalize(loaded.vectors(:,1:candidateCount), 1, 'norm');
    principalAnglesDeg(pointIndex,:) = loaded.principalAnglesRank4Deg;
end

trackedLambda = complex(zeros(pointCount, branchCount));
trackedCandidateIndex = zeros(pointCount, branchCount);
selectedMask = false(pointCount, branchCount);
trackedLambda(1,:) = lambdaByPoint{1}(1:branchCount);
trackedCandidateIndex(1,:) = 1:branchCount;
selectedMask(1,1:4) = true;
previousVectors = vectorsByPoint{1}(:,1:branchCount);

for pointIndex = 2:pointCount
    candidates = vectorsByPoint{pointIndex};
    candidateLambda = lambdaByPoint{pointIndex}(1:candidateCount);
    overlapCost = 1 - abs(previousVectors' * candidates);
    eigenvalueCost = abs(real(trackedLambda(pointIndex-1,:)).' - ...
        real(candidateLambda).');
    eigenvalueCost = min(eigenvalueCost / 0.02, 1);
    totalCost = 0.75 * overlapCost + 0.25 * eigenvalueCost;
    pairs = matchpairs(totalCost, 2);
    pairs = sortrows(pairs, 1);
    if size(pairs,1) ~= branchCount || any(pairs(:,1)' ~= 1:branchCount)
        error('Figure5D:BranchAssignment', ...
            'Could not assign all tracked branches at point %d.', pointIndex);
    end
    assigned = pairs(:,2)';
    trackedCandidateIndex(pointIndex,:) = assigned;
    trackedLambda(pointIndex,:) = candidateLambda(assigned);
    selectedMask(pointIndex,:) = assigned <= 4;
    previousVectors = candidates(:,assigned);
end

exchangeRows = find(any(selectedMask(2:end,:) ~= ...
    selectedMask(1:end-1,:), 2)) + 1;
tauMs = 10.3402405296839;
generatorBranchReal = (real(trackedLambda)-1) / tauMs;
insetData = struct('ArcLength',arcLength, ...
    'PrincipalAnglesDeg',principalAnglesDeg, ...
    'GeneratorBranchReal',generatorBranchReal, ...
    'SelectedMask',selectedMask, ...
    'TrackedCandidateIndex',trackedCandidateIndex, ...
    'ExchangeRows',exchangeRows);
save(fullfile(thisRoot, 'subspace_s0025', ...
    'figure_5D_inset_data.mat'), 'insetData');
fprintf('Tracked exchange rows: %s\n', mat2str(exchangeRows(:)'));
end
