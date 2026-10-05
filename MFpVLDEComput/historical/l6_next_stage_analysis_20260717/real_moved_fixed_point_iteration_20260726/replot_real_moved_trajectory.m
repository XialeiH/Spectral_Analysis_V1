function replot_real_moved_trajectory(setupFile, resultFile, summaryFile, outputFile)
% Replot the baseline and ten validated regular-iteration snapshots.

setupData = load(setupFile, 'setup');
resultData = load(resultFile, 'sampleIterations', 'sampleStates', ...
    'hcFromBaseline', 'hcToMoved');
summary = readtable(summaryFile, 'FileType', 'text', 'Delimiter', '\t');

context = setupData.setup.Context;
baseline = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
wC = context.CWeight;
snapshotIterations = resultData.sampleIterations(2:end);
snapshotStates = resultData.sampleStates(:,2:end);
hcBase = resultData.hcFromBaseline(2:end);
hcMoved = resultData.hcToMoved(2:end);

populationNames = {'S','C','I','E'};
stateCount = 1+numel(snapshotIterations);
allMaps = cell(4,stateCount);
allMaps(:,1) = local_state_maps(baseline,mapSize,wC);
for column = 2:stateCount
    allMaps(:,column) = local_state_maps( ...
        snapshotStates(:,column-1),mapSize,wC);
end
limits = local_population_limits(allMaps);

fig = figure('Visible','off','Color','w','Position',[20 20 4600 1750]);
layout = tiledlayout(fig,4,stateCount, ...
    'TileSpacing','compact','Padding','compact');
for row = 1:4
    for column = 1:stateCount
        ax = nexttile(layout,(row-1)*stateCount+column);
        imagesc(ax,allMaps{row,column});
        axis(ax,'image');
        axis(ax,'off');
        clim(ax,limits(row,:));
        if row == 1
            if column == 1
                title(ax,sprintf('Original baseline\nHC=0'), ...
                    'Interpreter','none','FontSize',13);
            else
                title(ax,sprintf(['iteration %d\nHC(base)=%.4g\n' ...
                    'HC(target)=%.3g'],snapshotIterations(column-1), ...
                    hcBase(column-1),hcMoved(column-1)),'FontSize',12);
            end
        end
        if column == 1
            text(ax,-0.08,0.5,populationNames{row},'Units','normalized', ...
                'HorizontalAlignment','right','FontWeight','bold','FontSize',18);
        elseif column == stateCount
            cb = colorbar(ax,'eastoutside');
            cb.FontSize = 13;
        end
    end
end
colormap(fig,parula(256));

titleText = sprintf(['%s regular moved-fixed-point iteration snapshots: ' ...
    'beta=%+.3f, gain=%.3f\nmax Re lambda=%.6f, p=%.2f, ' ...
    'final residual=%.2e'],char(summary.pathway(1)),summary.beta, ...
    summary.pathwayGain,summary.maxRealLambda,summary.relaxationP, ...
    summary.fixedPointResidual);
title(layout,titleText,'Interpreter','none','FontWeight','bold','FontSize',20);

exportgraphics(fig,outputFile,'ContentType','image','Resolution',180);
close(fig);
end

function maps = local_state_maps(state, mapSize, wC)
n = prod(mapSize);
maps = cell(4,1);
maps{1} = reshape(state(1:n),mapSize);
maps{2} = reshape(state(n+(1:n)),mapSize);
maps{3} = reshape(state(2*n+(1:n)),mapSize);
maps{4} = (1-wC)*maps{1}+wC*maps{2};
end

function limits = local_population_limits(allMaps)
limits = nan(4,2);
for row = 1:4
    values = [];
    for column = 1:size(allMaps,2)
        map = allMaps{row,column};
        values = [values; map(isfinite(map))]; %#ok<AGROW>
    end
    lower = min(values);
    upper = max(values);
    if upper <= lower
        padding = max(1e-6,abs(lower)*1e-6);
        lower = lower-padding;
        upper = upper+padding;
    end
    limits(row,:) = [lower upper];
end
end
