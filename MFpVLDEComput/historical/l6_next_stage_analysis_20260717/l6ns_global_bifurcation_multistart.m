function result = l6ns_global_bifurcation_multistart( ...
        setup,branchRoot,taskId,outputDir)
% Search for distinct fixed points at critical-side, w=0, and w=1 probes.

if taskId<1 || taskId>8
    error('taskId must be in 1:8.');
end
pathwayIndex = ceil(taskId/4);
probeIndex = mod(taskId-1,4)+1;
if pathwayIndex==1
    name = 'L6';
else
    name = 'Inhibition';
end
critical = setup.(name);
probeWeights = [critical.W-0.02 critical.W+0.02 0 1];
w = probeWeights(probeIndex);
context = setup.Context;
fixed = context.FixedPoint(:);
cfg = setup.Config;
if strcmp(name,'L6')
    phi = @(x,weight)l6ns_phi(x,weight,context,[1 1],1);
    persistentMatrix = l6ns_pathway_jacobian(setup.Pathway,1-w,1);
else
    phi = @(x,weight)l6ns_phi(x,0,context,[1 1],1-weight);
    persistentMatrix = l6ns_pathway_jacobian(setup.Pathway,1,1-w);
end

[branchSeeds,branchLabels] = local_branch_seeds(branchRoot,name,w,fixed,phi);
[seeds,labels] = local_seed_bank(fixed,critical,w,persistentMatrix,cfg, ...
    branchSeeds,branchLabels,taskId);

seedRows = {};
convergedStates = {};
for seedIndex = 1:numel(seeds)
    [state,diagnostic] = l6ns_fixed_point_newton(phi,fixed,w,seeds{seedIndex});
    seedRows(end+1,:) = {seedIndex,string(labels{seedIndex}),diagnostic.converged, ...
        diagnostic.iterations,diagnostic.residualNorm,diagnostic.relativeResidual, ...
        norm(state-fixed),min(state),max(state)}; %#ok<AGROW>
    if diagnostic.converged
        convergedStates{end+1,1} = state; %#ok<AGROW>
    end
    fprintf('%s w=%+.6f seed %d/%d %-24s residual %.3e converged=%d\n', ...
        name,w,seedIndex,numel(seeds),labels{seedIndex}, ...
        diagnostic.residualNorm,diagnostic.converged);
end
seedTable = cell2table(seedRows,'VariableNames',{'seedIndex','seedType', ...
    'converged','newtonIterations','residualNorm','relativeResidual', ...
    'distanceFromPersistent','minimumState','maximumState'});

[uniqueStates,assignments] = local_deduplicate(convergedStates,fixed);
rootRows = {};
for rootIndex = 1:numel(uniqueStates)
    state = uniqueStates{rootIndex};
    modes = l6ns_numeric_leading_modes(@(x)phi(x,w),state,4,cfg);
    residualNorm = norm(phi(state,w)-state);
    rootRows(end+1,:) = {rootIndex,norm(state-fixed),residualNorm, ...
        residualNorm/max(1,norm(fixed)),max(real(modes.Lambda)), ...
        min(state),max(state),sum(assignments==rootIndex)}; %#ok<AGROW>
end
rootTable = cell2table(rootRows,'VariableNames',{'rootIndex', ...
    'distanceFromPersistent','residualNorm','relativeResidual','maxRealLambda', ...
    'minimumState','maximumState','attractionCount'});

if ~exist(outputDir,'dir'); mkdir(outputDir); end
stem = sprintf('%s_probe_%d_w_%+.6f',lower(name),probeIndex,w);
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
stem = strrep(stem,'.','p');
writetable(seedTable,fullfile(outputDir,[stem '_seeds.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(rootTable,fullfile(outputDir,[stem '_roots.tsv']), ...
    'FileType','text','Delimiter','\t');
result = struct('Name',name,'ProbeIndex',probeIndex,'FreezeWeight',w, ...
    'SeedTable',seedTable,'RootTable',rootTable,'RootStates',{uniqueStates}, ...
    'ConvergedAssignments',assignments);
save(fullfile(outputDir,[stem '_result.mat']),'result','-v7.3');
end

function [states,labels] = local_branch_seeds(branchRoot,name,w,fixed,phi)
states = {};
labels = {};
for branchSign = [-1 1]
    stem = sprintf('%s_sign_%+d_result.mat',lower(name),branchSign);
    stem = strrep(stem,'+','p');
    stem = strrep(stem,'-','m');
    fileName = fullfile(branchRoot,stem);
    if ~isfile(fileName); continue; end
    loaded = load(fileName,'result');
    tableOut = loaded.result.Table;
    valid = find(tableOut.converged);
    bracket = [];
    for index = 1:numel(valid)-1
        first = valid(index);
        second = valid(index+1);
        if (tableOut.freezeWeight(first)-w)*(tableOut.freezeWeight(second)-w)<=0
            bracket = [first second]; %#ok<AGROW>
            break
        end
    end
    if isempty(bracket); continue; end
    w1 = tableOut.freezeWeight(bracket(1));
    w2 = tableOut.freezeWeight(bracket(2));
    fraction = (w-w1)/(w2-w1);
    guess = loaded.result.States{bracket(1)}+fraction* ...
        (loaded.result.States{bracket(2)}-loaded.result.States{bracket(1)});
    [state,diagnostic] = l6ns_fixed_point_newton(phi,fixed,w,guess);
    if diagnostic.converged
        states{end+1,1} = state; %#ok<AGROW>
        labels{end+1,1} = sprintf('continued_%+d',branchSign); %#ok<AGROW>
    end
end
end

function [seeds,labels] = local_seed_bank(fixed,critical,w,matrix,cfg, ...
        branchSeeds,branchLabels,seed)
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
seeds = [{fixed};branchSeeds(:)];
labels = [{'persistent'};branchLabels(:)];

normalAmplitude = sqrt(max(0,-critical.Beta*(w-critical.W)/critical.Cubic));
if ~isfinite(normalAmplitude) || normalAmplitude<stateScale*0.1
    normalAmplitude = stateScale*4;
end
for branchSign = [-1 1]
    for factor = [0.25 0.6 1 1.4]
        seeds{end+1,1} = fixed+branchSign*factor*normalAmplitude*critical.Right; %#ok<AGROW>
        labels{end+1,1} = sprintf('critical_%+d_x%.2g',branchSign,factor); %#ok<AGROW>
    end
end

modes = l6ns_eigenpairs(matrix,8,'largestreal',cfg);
basis = real(modes.Right);
[basis,~] = qr(basis,0);
rng(cfg.RandomSeed+seed);
radius = normalAmplitude;
if ~isempty(branchSeeds)
    distances = cellfun(@(x)norm(x-fixed),branchSeeds);
    radius = max(radius,mean(distances));
end
for randomIndex = 1:12
    coefficients = randn(size(basis,2),1);
    direction = basis*coefficients;
    direction = direction/max(norm(direction),eps);
    factor = [0.35 0.7 1.1 1.6];
    factor = factor(mod(randomIndex-1,numel(factor))+1);
    seeds{end+1,1} = fixed+factor*radius*direction; %#ok<AGROW>
    labels{end+1,1} = sprintf('leading_subspace_%02d',randomIndex); %#ok<AGROW>
end
seeds{end+1,1} = 0.5*fixed;
labels{end+1,1} = 'half_persistent';
seeds{end+1,1} = 1.5*fixed;
labels{end+1,1} = 'onehalf_persistent';
end

function [uniqueStates,assignments] = local_deduplicate(states,fixed)
uniqueStates = {};
assignments = zeros(numel(states),1);
tolerance = 2e-6*max(1,norm(fixed));
for stateIndex = 1:numel(states)
    state = states{stateIndex};
    match = 0;
    for rootIndex = 1:numel(uniqueStates)
        if norm(state-uniqueStates{rootIndex})<tolerance
            match = rootIndex;
            break
        end
    end
    if match==0
        uniqueStates{end+1,1} = state; %#ok<AGROW>
        match = numel(uniqueStates);
    end
    assignments(stateIndex) = match;
end
end
