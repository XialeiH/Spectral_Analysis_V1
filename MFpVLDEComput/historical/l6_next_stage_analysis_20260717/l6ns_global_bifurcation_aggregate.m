function result = l6ns_global_bifurcation_aggregate(setup,branchRoot,multistartRoot,outputDir)
% Aggregate wide continuation, stability samples, and the multistart census.

if ~exist(outputDir,'dir'); mkdir(outputDir); end
branch = struct();
for name = {'L6','Inhibition'}
    field = name{1};
    branch.(field) = cell(2,1);
    for signIndex = 1:2
        branchSign = [-1 1];
        stem = sprintf('%s_sign_%+d_result.mat',lower(field),branchSign(signIndex));
        stem = strrep(stem,'+','p');
        stem = strrep(stem,'-','m');
        loaded = load(fullfile(branchRoot,stem),'result');
        branch.(field){signIndex} = loaded.result;
    end
end

multistartFiles = dir(fullfile(multistartRoot,'*_result.mat'));
if numel(multistartFiles)~=8
    error('Expected 8 multistart results, found %d.',numel(multistartFiles));
end
multistart = cell(numel(multistartFiles),1);
for index = 1:numel(multistartFiles)
    loaded = load(fullfile(multistartFiles(index).folder,multistartFiles(index).name),'result');
    multistart{index} = loaded.result;
end

stabilitySamples = local_stability_samples(setup,branch);
rootCensus = local_root_census(setup,multistart,branch);
symmetryPairing = local_symmetry_pairing(setup,multistart);
persistentStability = local_persistent_stability(setup);

writetable(stabilitySamples,fullfile(outputDir,'secondary_branch_stability_samples.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(rootCensus,fullfile(outputDir,'fixed_point_root_census.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(symmetryPairing,fullfile(outputDir,'secondary_branch_symmetry_pairing.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(persistentStability,fullfile(outputDir,'persistent_branch_stability_wide.tsv'), ...
    'FileType','text','Delimiter','\t');
copyfile(fullfile(setup.OutputDir,'spatial_symmetry_audit.tsv'), ...
    fullfile(outputDir,'spatial_symmetry_audit.tsv'));

local_plot_pathway(setup,branch.L6,multistart,stabilitySamples, ...
    persistentStability,'L6',outputDir);
local_plot_pathway(setup,branch.Inhibition,multistart,stabilitySamples, ...
    persistentStability,'Inhibition',outputDir);
local_plot_root_counts(rootCensus,outputDir);
local_plot_symmetry(setup.SymmetryAudit,symmetryPairing,outputDir);
extraBranchFile = fullfile(fileparts(branchRoot),'extra_branches', ...
    'inhibition_nearsilent_branch_result.mat');
if isfile(extraBranchFile)
    extraLoaded = load(extraBranchFile,'result');
    local_plot_inhibition_combined(setup,branch.Inhibition,stabilitySamples, ...
        persistentStability,extraLoaded.result,outputDir);
end

result = struct('SetupSummary',local_setup_summary(setup), ...
    'Branch',branch,'Multistart',{multistart}, ...
    'SecondaryStability',stabilitySamples,'RootCensus',rootCensus, ...
    'SymmetryPairing',symmetryPairing,'PersistentStability',persistentStability);
save(fullfile(outputDir,'global_bifurcation_extension_result.mat'),'result','-v7.3');
end

function tableOut = local_stability_samples(setup,branch)
rows = {};
for name = {'L6','Inhibition'}
    field = name{1};
    context = setup.Context;
    if strcmp(field,'L6')
        phi = @(x,w)l6ns_phi(x,w,context,[1 1],1);
    else
        phi = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
    end
    for signIndex = 1:2
        item = branch.(field){signIndex};
        valid = find(item.Table.converged);
        samplePositions = unique(round(linspace(1,numel(valid),min(9,numel(valid)))));
        for position = samplePositions
            rowIndex = valid(position);
            state = item.States{rowIndex};
            w = item.Table.freezeWeight(rowIndex);
            modes = l6ns_numeric_leading_modes(@(x)phi(x,w),state,4,setup.Config);
            rows(end+1,:) = {string(field),item.Sign,w, ...
                item.Table.normalizedAmplitude(rowIndex),norm(state-context.FixedPoint), ...
                max(real(modes.Lambda)),min(state),max(state)}; %#ok<AGROW>
        end
    end
end
tableOut = cell2table(rows,'VariableNames',{'pathway','branchSign','freezeWeight', ...
    'normalizedAmplitude','branchDistance','maxRealLambda','minimumState','maximumState'});
end

function tableOut = local_root_census(setup,multistart,branch)
rows = {};
fixed = setup.Context.FixedPoint(:);
dedupTolerance = 2e-5*max(1,norm(fixed));
for resultIndex = 1:numel(multistart)
    item = multistart{resultIndex};
    critical = setup.(item.Name);
    states = {};
    sources = strings(0,1);
    attractionCounts = zeros(0,1);
    if strcmp(item.Name,'L6')
        phi = @(x,w)l6ns_phi(x,w,setup.Context,[1 1],1);
    else
        phi = @(x,w)l6ns_phi(x,0,setup.Context,[1 1],1-w);
    end
    for rawIndex = 1:numel(item.RootStates)
        [candidate,polish] = l6ns_fixed_point_polish(phi,fixed, ...
            item.FreezeWeight,item.RootStates{rawIndex});
        if ~polish.converged
            candidate = item.RootStates{rawIndex};
            source = "multistart_unpolished";
        else
            source = "multistart";
        end
        [states,sources,attractionCounts] = local_add_root(states,sources, ...
            attractionCounts,candidate,source, ...
            item.RootTable.attractionCount(rawIndex),dedupTolerance);
    end
    for signIndex = 1:2
        branchItem = branch.(item.Name){signIndex};
        valid = find(branchItem.Table.converged);
        for segmentIndex = 1:numel(valid)-1
            first = valid(segmentIndex);
            second = valid(segmentIndex+1);
            firstW = branchItem.Table.freezeWeight(first);
            secondW = branchItem.Table.freezeWeight(second);
            if (firstW-item.FreezeWeight)*(secondW-item.FreezeWeight)>0 || firstW==secondW
                continue
            end
            fraction = (item.FreezeWeight-firstW)/(secondW-firstW);
            guess = branchItem.States{first}+fraction* ...
                (branchItem.States{second}-branchItem.States{first});
            [candidate,diagnostic] = l6ns_fixed_point_polish( ...
                phi,fixed,item.FreezeWeight,guess);
            if ~diagnostic.converged
                continue
            end
            [states,sources,attractionCounts] = local_add_root(states,sources, ...
                attractionCounts,candidate,"continuation",0,dedupTolerance);
        end
    end

    seedCount = height(item.SeedTable);
    convergedSeedCount = sum(item.SeedTable.converged);
    for rowIndex = 1:numel(states)
        state = states{rowIndex};
        residualNorm = norm(phi(state,item.FreezeWeight)-state);
        if sources(rowIndex)=="continuation"
            modes = l6ns_numeric_leading_modes(@(x)phi(x,item.FreezeWeight), ...
                state,4,setup.Config);
            maxRealLambda = max(real(modes.Lambda));
        else
            rawDistances = cellfun(@(x)norm(state-x),item.RootStates);
            [~,nearestRaw] = min(rawDistances);
            maxRealLambda = item.RootTable.maxRealLambda(nearestRaw);
        end
        attractionCount = attractionCounts(rowIndex);
        signedAmplitude = real(critical.Left'*(state-fixed));
        rows(end+1,:) = {string(item.Name),item.ProbeIndex,item.FreezeWeight, ... %#ok<AGROW>
            rowIndex,sources(rowIndex),signedAmplitude,norm(state-fixed), ...
            residualNorm,residualNorm/max(1,norm(fixed)),maxRealLambda, ...
            maxRealLambda<1,min(state),max(state),attractionCount,seedCount, ...
            convergedSeedCount};
    end
end
tableOut = cell2table(rows,'VariableNames',{'pathway','probeIndex','freezeWeight', ...
    'rootIndex','source','signedCriticalAmplitude','distanceFromPersistent','residualNorm', ...
    'relativeResidual','maxRealLambda','stable','minimumState','maximumState', ...
    'seedAttractionCount','seedCount','convergedSeedCount'});
tableOut = sortrows(tableOut,{'pathway','probeIndex','signedCriticalAmplitude'});
end

function [states,sources,counts] = local_add_root( ...
        states,sources,counts,candidate,source,count,tolerance)
match = 0;
for rootIndex = 1:numel(states)
    if norm(candidate-states{rootIndex})<tolerance
        match = rootIndex;
        break
    end
end
if match==0
    states{end+1,1} = candidate;
    sources(end+1,1) = source;
    counts(end+1,1) = count;
else
    counts(match) = counts(match)+count;
    if sources(match)~="multistart" && source=="multistart"
        sources(match) = source;
    end
end
end

function tableOut = local_symmetry_pairing(setup,multistart)
rows = {};
fixed = setup.Context.FixedPoint(:);
mapSize = setup.Context.MapSize;
for name = {'L6','Inhibition'}
    field = name{1};
    auditName = field;
    if strcmp(field,'Inhibition'); auditName = 'I'; end
    audit = setup.SymmetryAudit(setup.SymmetryAudit.pathway==string(auditName),:);
    score = audit.fixedPointInvarianceError+audit.criticalModeOddError+ ...
        audit.mapEquivarianceError;
    [~,best] = min(score);
    rowShift = audit.rowShift(best);
    columnShift = audit.columnShift(best);
    matching = multistart(cellfun(@(x)strcmp(x.Name,field),multistart));
    for resultIndex = 1:numel(matching)
        item = matching{resultIndex};
        amplitudes = cellfun(@(x)real(setup.(field).Left'*(x-fixed)),item.RootStates);
        distances = item.RootTable.distanceFromPersistent;
        secondary = find(distances>2e-5*max(1,norm(fixed)));
        if isempty(secondary)
            rows(end+1,:) = {string(field),item.ProbeIndex,item.FreezeWeight, ... %#ok<AGROW>
                rowShift,columnShift,0,0,NaN,NaN,NaN,NaN};
            continue
        end
        for source = secondary(:)'
            shifted = local_shift_state(item.RootStates{source},mapSize,rowShift,columnShift);
            relativeErrors = inf(size(secondary));
            for candidateIndex = 1:numel(secondary)
                candidate = secondary(candidateIndex);
                relativeErrors(candidateIndex) = norm(shifted-item.RootStates{candidate}) / ...
                    max(norm(item.RootStates{source}-fixed),eps);
            end
            [pairingError,bestCandidate] = min(relativeErrors);
            target = secondary(bestCandidate);
            distanceRatio = distances(source)/max(distances(target),eps);
            rows(end+1,:) = {string(field),item.ProbeIndex,item.FreezeWeight, ... %#ok<AGROW>
                rowShift,columnShift,source,target,amplitudes(source),amplitudes(target), ...
                pairingError,distanceRatio};
        end
    end
end
tableOut = cell2table(rows,'VariableNames',{'pathway','probeIndex','freezeWeight', ...
    'rowShift','columnShift','sourceRootIndex','matchedRootIndex', ...
    'sourceSignedAmplitude','matchedSignedAmplitude', ...
    'symmetryPairingRelativeError','branchDistanceRatio'});
tableOut = sortrows(tableOut,{'pathway','probeIndex','sourceRootIndex'});
end

function tableOut = local_persistent_stability(setup)
wGrid = unique([-0.30:0.025:1.10,setup.L6.W,setup.Inhibition.W]);
rows = {};
for name = {'L6','Inhibition'}
    field = name{1};
    for w = wGrid
        if strcmp(field,'L6')
            matrix = l6ns_pathway_jacobian(setup.Pathway,1-w,1);
        else
            matrix = l6ns_pathway_jacobian(setup.Pathway,1,1-w);
        end
        alpha = l6ns_max_real(matrix,setup.Config);
        rows(end+1,:) = {string(field),w,alpha,alpha<1}; %#ok<AGROW>
    end
end
tableOut = cell2table(rows,'VariableNames', ...
    {'pathway','freezeWeight','maxRealLambda','stable'});
end

function summary = local_setup_summary(setup)
summary = table(["L6";"Inhibition"],[setup.L6.W;setup.Inhibition.W], ...
    [setup.L6.Beta;setup.Inhibition.Beta], ...
    [setup.L6.Quadratic;setup.Inhibition.Quadratic], ...
    [setup.L6.Cubic;setup.Inhibition.Cubic], ...
    [setup.L6.SpectralSeparation;setup.Inhibition.SpectralSeparation], ...
    'VariableNames',{'pathway','criticalFreezeWeight','beta', ...
    'quadraticCoefficient','cubicCoefficient','spectralSeparation'});
end

function local_plot_pathway(setup,branches,multistart,stability,persistentTable,name,outputDir)
critical = setup.(name);
items = multistart(cellfun(@(x)strcmp(x.Name,name),multistart));
fig = figure('Visible','off','Color','w','Position',[80 60 1180 900]);
layout = tiledlayout(2,1,'TileSpacing','loose','Padding','loose');

ax1 = nexttile(layout); hold(ax1,'on');
persistentRows = persistentTable.pathway==string(name);
local_plot_stability_segments(ax1,persistentTable.freezeWeight(persistentRows), ...
    zeros(sum(persistentRows),1),persistentTable.stable(persistentRows),'k',2.0);
colors = [0.10 0.45 0.80;0.85 0.33 0.10];
for signIndex = 1:2
    tableOut = branches{signIndex}.Table;
    valid = tableOut.converged;
    plot(ax1,tableOut.freezeWeight(valid),tableOut.normalizedAmplitude(valid), ...
        '-','Color',colors(signIndex,:),'LineWidth',1.7, ...
        'DisplayName',sprintf('secondary sign %+d',branches{signIndex}.Sign));
end
for itemIndex = 1:numel(items)
    item = items{itemIndex};
    amplitudes = cellfun(@(x)real(critical.Left'*(x-setup.Context.FixedPoint)), ...
        item.RootStates) / max(1,norm(setup.Context.FixedPoint)/sqrt(numel(setup.Context.FixedPoint)));
    stable = item.RootTable.maxRealLambda<1;
    scatter(ax1,repmat(item.FreezeWeight,size(amplitudes)),amplitudes,52, ...
        stable,'filled','MarkerEdgeColor',[0.15 0.15 0.15], ...
        'HandleVisibility','off');
end
xline(ax1,critical.W,'b:','LineWidth',1.3,'DisplayName','critical weight');
xline(ax1,0,':','Color',[0.4 0.4 0.4],'HandleVisibility','off');
xline(ax1,1,':','Color',[0.4 0.4 0.4],'HandleVisibility','off');
grid(ax1,'on'); xlabel(ax1,'freeze weight w');
ylabel(ax1,'signed critical amplitude / state scale');
title(ax1,sprintf('%s: wide fixed-point continuation and root census',name));
legend(ax1,'Location','best'); colorbar(ax1,'Ticks',[0 1], ...
    'TickLabels',{'unstable','stable'});
clim(ax1,[0 1]);

ax2 = nexttile(layout); hold(ax2,'on');
plot(ax2,persistentTable.freezeWeight(persistentRows), ...
    persistentTable.maxRealLambda(persistentRows),'k-','LineWidth',1.8, ...
    'DisplayName','persistent branch');
sampleRows = stability.pathway==string(name);
for signValue = [-1 1]
    rows = sampleRows & stability.branchSign==signValue;
    color = colors((signValue+3)/2,:);
    plot(ax2,stability.freezeWeight(rows),stability.maxRealLambda(rows),'o-', ...
        'Color',color,'MarkerFaceColor',color,'LineWidth',1.3, ...
        'DisplayName',sprintf('secondary sign %+d',signValue));
end
yline(ax2,1,'r--','LineWidth',1.2,'DisplayName','stability boundary');
xline(ax2,critical.W,'b:','LineWidth',1.2,'HandleVisibility','off');
grid(ax2,'on'); xlabel(ax2,'freeze weight w'); ylabel(ax2,'max Re \lambda(D\Phi)');
title(ax2,'Persistent stability and sampled secondary-branch stability');
legend(ax2,'Location','best');
set([ax1 ax2],'FontSize',11);
stem = sprintf('%s_global_bifurcation_extension',lower(name));
l6ns_save_figure(fig,outputDir,stem); close(fig);
end

function local_plot_stability_segments(ax,x,y,stable,color,width)
[x,order] = sort(x);
y = y(order);
stable = stable(order);
for index = 1:numel(x)-1
    if stable(index) && stable(index+1); style='-'; else; style='--'; end
    plot(ax,x(index:index+1),y(index:index+1),style,'Color',color, ...
        'LineWidth',width,'HandleVisibility','off');
end
plot(ax,NaN,NaN,'k-','LineWidth',width,'DisplayName','persistent stable');
plot(ax,NaN,NaN,'k--','LineWidth',width,'DisplayName','persistent unstable');
end

function local_plot_root_counts(census,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1050 430]);
layout = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');
for name = {'L6','Inhibition'}
    ax = nexttile(layout);
    rows = census.pathway==string(name{1});
    weights = unique(census.freezeWeight(rows),'stable');
    counts = arrayfun(@(w)sum(rows & census.freezeWeight==w),weights);
    stableCounts = arrayfun(@(w)sum(rows & census.freezeWeight==w & census.stable),weights);
    unstableCounts = counts-stableCounts;
    bar(ax,weights,[stableCounts unstableCounts],'stacked'); grid(ax,'on');
    xlabel(ax,'freeze weight w'); ylabel(ax,'distinct converged roots');
    title(ax,sprintf('%s polished root census (multistart + continuation)',name{1}));
    legend(ax,{'stable','unstable'},'Location','best');
end
l6ns_save_figure(fig,outputDir,'global_bifurcation_root_census'); close(fig);
end

function local_plot_symmetry(audit,pairing,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1100 430]);
layout = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');
for name = {'L6','Inhibition'}
    ax = nexttile(layout); hold(ax,'on');
    auditName = name{1};
    if strcmp(auditName,'Inhibition'); auditName = 'I'; end
    rows = audit.pathway==string(auditName);
    values = [audit.fixedPointInvarianceError(rows),audit.criticalModeOddError(rows), ...
        audit.mapEquivarianceError(rows)];
    score = sum(values,2);
    [~,order] = sort(score);
    order = order(1:min(6,numel(order)));
    semilogy(ax,1:numel(order),max(values(order,:),1e-18),'-o','LineWidth',1.2);
    pairRows = pairing.pathway==string(name{1}) & ...
        isfinite(pairing.symmetryPairingRelativeError);
    if any(pairRows)
        yline(ax,median(pairing.symmetryPairingRelativeError(pairRows)), ...
            'k--','LineWidth',1.2,'Label','median branch-pair error');
    end
    grid(ax,'on'); xlabel(ax,'best symmetry candidates'); ylabel(ax,'relative error');
    title(ax,sprintf('%s spatial-symmetry audit',name{1}));
    legend(ax,{'fixed point','critical mode oddness','map equivariance'}, ...
        'Location','best');
end
l6ns_save_figure(fig,outputDir,'global_bifurcation_symmetry_audit'); close(fig);
end

function local_plot_inhibition_combined(setup,branches,stability,persistentTable, ...
        extraBranch,outputDir)
fixedNorm = norm(setup.Context.FixedPoint(:));
fig = figure('Visible','off','Color','w','Position',[80 60 1180 900]);
layout = tiledlayout(2,1,'TileSpacing','loose','Padding','loose');
colors = [0.10 0.45 0.80;0.85 0.33 0.10];

ax1 = nexttile(layout); hold(ax1,'on');
plot(ax1,[-0.3 1.1],[0 0],'k-','LineWidth',1.8,'DisplayName','persistent branch');
for signIndex = 1:2
    tableOut = branches{signIndex}.Table;
    valid = tableOut.converged;
    plot(ax1,tableOut.freezeWeight(valid),tableOut.branchDistance(valid)/fixedNorm, ...
        '-','Color',colors(signIndex,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('pitchfork family sign %+d',branches{signIndex}.Sign));
end
plot(ax1,extraBranch.Table.freezeWeight,extraBranch.Table.relativeDistance, ...
    'm-','LineWidth',2.0,'DisplayName','near-silent saddle-node family');
xline(ax1,setup.Inhibition.W,'r:','LineWidth',1.2,'DisplayName','pitchfork critical weight');
grid(ax1,'on'); xlabel(ax1,'inhibition freeze weight w_I');
ylabel(ax1,'||f-f^*||_2 / ||f^*||_2');
title(ax1,'Inhibition: persistent, pitchfork, and near-silent fixed-point families');
legend(ax1,'Location','best');

ax2 = nexttile(layout); hold(ax2,'on');
rows = persistentTable.pathway=="Inhibition";
plot(ax2,persistentTable.freezeWeight(rows),persistentTable.maxRealLambda(rows), ...
    'k-','LineWidth',1.8,'DisplayName','persistent branch');
sampleRows = stability.pathway=="Inhibition";
for signValue = [-1 1]
    rows = sampleRows & stability.branchSign==signValue;
    plot(ax2,stability.freezeWeight(rows),stability.maxRealLambda(rows),'o-', ...
        'Color',colors((signValue+3)/2,:),'MarkerFaceColor',colors((signValue+3)/2,:), ...
        'LineWidth',1.2,'DisplayName',sprintf('pitchfork sign %+d',signValue));
end
plot(ax2,extraBranch.Stability.freezeWeight,extraBranch.Stability.maxRealLambda, ...
    'mo-','MarkerFaceColor','m','LineWidth',1.5,'DisplayName','near-silent family');
yline(ax2,1,'r--','LineWidth',1.2,'DisplayName','stability boundary');
xline(ax2,setup.Inhibition.W,'r:','LineWidth',1.2,'HandleVisibility','off');
grid(ax2,'on'); xlabel(ax2,'inhibition freeze weight w_I');
ylabel(ax2,'max Re \lambda(D\Phi)'); title(ax2,'Sampled stability of all inhibition families');
legend(ax2,'Location','best'); set([ax1 ax2],'FontSize',11);
l6ns_save_figure(fig,outputDir,'inhibition_all_fixed_point_families'); close(fig);
end

function shifted = local_shift_state(state,mapSize,rowShift,columnShift)
n = prod(mapSize);
shifted = zeros(size(state));
for population = 1:3
    indices = (population-1)*n+(1:n);
    map = reshape(state(indices),mapSize);
    shifted(indices) = reshape(circshift(map,[rowShift columnShift]),[],1);
end
end
