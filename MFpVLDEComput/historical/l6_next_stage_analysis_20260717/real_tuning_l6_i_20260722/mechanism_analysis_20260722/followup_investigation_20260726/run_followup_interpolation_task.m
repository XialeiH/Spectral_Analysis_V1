function run_followup_interpolation_task()
% Test causal state interpolation, population factorials, and mean/pattern shifts.

paths = followup_initialize();
loaded = load(paths.SetupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
specs = followup_specs('interpolation');
taskId = local_task_id();
spec = specs(taskId,:);
[gain6,gainI] = local_gains(spec.pathway,spec.beta);
moved = followup_branch_state(paths,spec.pathway,spec.beta,"low",context);
n = numel(baseline)/3;

alphas = (0:0.1:1)';
interpolationCells = cell(numel(alphas),1);
slopeCells = cell(numel(alphas)*9,1);
slopeIndex = 0;
for index = 1:numel(alphas)
    alpha = alphas(index);
    state = (1-alpha)*baseline+alpha*moved;
    [J,~,details] = mechanism_jacobian_components( ...
        state,context,gain6,gainI);
    maxReal = mechanism_max_real(J);
    interpolationCells{index} = table(taskId,spec.pathway,spec.beta,alpha, ...
        maxReal,norm(state-baseline)/norm(baseline),min(state),max(state), ...
        'VariableNames',{'taskId','pathway','beta','alpha','maxRealLambda', ...
        'relativeStateShift','minimumRateHz','maximumRateHz'});
    localRows = local_slope_rows(taskId,spec.pathway,spec.beta,alpha,details);
    for row = 1:height(localRows)
        slopeIndex = slopeIndex+1;
        slopeCells{slopeIndex} = localRows(row,:);
    end
end
interpolationTable = vertcat(interpolationCells{:});
slopeTable = vertcat(slopeCells{1:slopeIndex});

factorialCells = cell(8,1);
factorialValues = nan(8,1);
for mask = 0:7
    state = baseline;
    movedFlags = false(1,3);
    for population = 1:3
        movedFlags(population) = bitget(mask,population)>0;
        if movedFlags(population)
            indices = (population-1)*n+(1:n);
            state(indices) = moved(indices);
        end
    end
    maxReal = mechanism_max_real(real_tuning_true_jacobian( ...
        state,context,gain6,gainI));
    factorialValues(mask+1) = maxReal;
    factorialCells{mask+1} = table(taskId,spec.pathway,spec.beta,mask, ...
        movedFlags(1),movedFlags(2),movedFlags(3),maxReal, ...
        'VariableNames',{'taskId','pathway','beta','subsetMask', ...
        'movedS','movedC','movedI','maxRealLambda'});
end
factorialTable = vertcat(factorialCells{:});
shapley = local_shapley(factorialValues);
shapleyTable = table(repmat(taskId,3,1),repmat(spec.pathway,3,1), ...
    repmat(spec.beta,3,1),["S";"C";"I"],shapley(:), ...
    'VariableNames',{'taskId','pathway','beta','population','shapleyContribution'});

delta = moved-baseline;
deltaMean = zeros(size(delta));
for population = 1:3
    indices = (population-1)*n+(1:n);
    deltaMean(indices) = mean(delta(indices));
end
deltaPattern = delta-deltaMean;
states = {baseline,baseline+deltaMean,baseline+deltaPattern,moved};
labels = ["baseline";"mean only";"pattern only";"full moved"];
meanPatternCells = cell(4,1);
for index = 1:4
    maxReal = mechanism_max_real(real_tuning_true_jacobian( ...
        states{index},context,gain6,gainI));
    meanPatternCells{index} = table(taskId,spec.pathway,spec.beta,labels(index), ...
        maxReal,norm(states{index}-baseline)/norm(baseline), ...
        'VariableNames',{'taskId','pathway','beta','stateComponent', ...
        'maxRealLambda','relativeStateShift'});
end
meanPatternTable = vertcat(meanPatternCells{:});

outDir = fullfile(paths.FollowupOutput,'interpolation_tasks');
if ~exist(outDir,'dir'); mkdir(outDir); end
stem = sprintf('%03d',taskId);
writetable(interpolationTable,fullfile(outDir,['interpolation_' stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(slopeTable,fullfile(outDir,['interpolation_slopes_' stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(factorialTable,fullfile(outDir,['factorial_' stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(shapleyTable,fullfile(outDir,['shapley_' stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(meanPatternTable,fullfile(outDir,['mean_pattern_' stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outDir,['interpolation_' stem '.mat']), ...
    'interpolationTable','slopeTable','factorialTable','shapleyTable', ...
    'meanPatternTable','baseline','moved','-v7.3');
fprintf('%s beta=%+.6f interpolation rho %.6f -> %.6f; Shapley S/C/I %+.5g %+.5g %+.5g\n', ...
    spec.pathway,spec.beta,interpolationTable.maxRealLambda(1), ...
    interpolationTable.maxRealLambda(end),shapley(1),shapley(2),shapley(3));
end

function contributions = local_shapley(values)
contributions = zeros(3,1);
for population = 1:3
    for mask = 0:7
        if bitget(mask,population); continue; end
        subsetSize = sum(bitget(mask,1:3));
        weight = factorial(subsetSize)*factorial(2-subsetSize)/factorial(3);
        withPopulation = bitset(mask,population,1);
        contributions(population) = contributions(population)+weight* ...
            (values(withPopulation+1)-values(mask+1));
    end
end
end

function rows = local_slope_rows(taskId,pathway,beta,alpha,details)
n = numel(details.L6);
targets = ["S";"C";"I"];
channels = ["L4E";"L4I";"L6"];
valuesByChannel = {details.L4EGradient,details.L4IGradient,details.L6Gradient};
cells = cell(9,1); rowIndex=0;
for target=1:3
    indices=(target-1)*n+(1:n);
    for channel=1:3
        rowIndex=rowIndex+1;
        values=valuesByChannel{channel}(indices);
        q=quantile(values,[0.05 0.5 0.95]);
        cells{rowIndex}=table(taskId,pathway,beta,alpha,targets(target), ...
            channels(channel),q(1),q(2),q(3),max(values), ...
            mean(abs(values)<1e-6), ...
            'VariableNames',{'taskId','pathway','beta','alpha','target', ...
            'channel','q05','median','q95','maximum','fractionAbsBelow1eMinus6'});
    end
end
rows=vertcat(cells{:});
end

function taskId = local_task_id()
taskId=str2double(getenv('FOLLOWUP_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('SLURM_ARRAY_TASK_ID')); end
if ~isfinite(taskId); taskId=1; end
taskId=round(taskId);
end

function [gain6,gainI] = local_gains(pathway,beta)
gain6=1; gainI=1;
if pathway=="L6"; gain6=1+beta; else; gainI=1+beta; end
end
