function replot_real_moved_pair(setupFile, resultFile, summaryFile, outputFile)
% Replot a validated baseline-versus-moved fixed-point comparison.

setupData = load(setupFile, 'setup');
resultData = load(resultFile, 'moved');
summary = readtable(summaryFile, 'FileType', 'text', 'Delimiter', '\t');

context = setupData.setup.Context;
baseline = context.FixedPoint(:);
moved = resultData.moved(:);
mapSize = double(context.MapSize(:).');
wC = context.CWeight;

populationNames = {'S','C','I','E'};
states = {baseline, moved};
allMaps = cell(4,2);
for column = 1:2
    allMaps(:,column) = local_state_maps(states{column}, mapSize, wC);
end
limits = local_population_limits(allMaps);

fig = figure('Visible','off','Color','w','Position',[40 40 1150 1700]);
layout = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
for row = 1:4
    for column = 1:2
        ax = nexttile(layout,(row-1)*2+column);
        imagesc(ax,allMaps{row,column});
        axis(ax,'image');
        axis(ax,'off');
        clim(ax,limits(row,:));
        if row == 1
            if column == 1
                title(ax,'Original baseline fixed point, HC=0','FontSize',12);
            else
                title(ax,sprintf(['Moved fixed point after %d iterations\n' ...
                    'HC(base)=%.5g, residual=%.2e'], ...
                    summary.fixedPointIterations, ...
                    summary.movedHCnormFromOriginalBaseline, ...
                    summary.fixedPointResidual),'FontSize',12);
            end
        end
        if column == 1
            text(ax,-0.08,0.5,populationNames{row},'Units','normalized', ...
                'HorizontalAlignment','right','FontWeight','bold','FontSize',13);
        else
            cb = colorbar(ax,'eastoutside');
            cb.FontSize = 11;
        end
    end
end
colormap(fig,parula(256));

titleText = sprintf(['%s real moved-fixed-point iteration: beta=%+.3f, ' ...
    'gain=%.3f\nmax Re lambda=%.6f, p=%.2f'], ...
    char(summary.pathway(1)),summary.beta,summary.pathwayGain, ...
    summary.maxRealLambda,summary.relaxationP);
title(layout,titleText,'Interpreter','none','FontWeight','bold','FontSize',13);

exportgraphics(fig,outputFile,'ContentType','image','Resolution',220);
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
