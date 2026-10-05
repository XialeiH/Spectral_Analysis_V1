function make_Final_Figure_4D()
% Replot the selected HC=1 perturbation snapshots for report Figure 4E.

artifactRoot = fileparts(mfilename('fullpath'));
projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
sourceFile = fullfile(projectRoot, ...
    'Spectral_Analysis/matlab-inserting_into_CG_model', ...
    'l6_next_stage_analysis_20260717/perturbation_7_20260807', ...
    'all_figures_finite_time_sources/7.4_source_data.mat');
outputPdf = fullfile(projectRoot, ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 4D.refresh.pdf');
finalPdf = fullfile(projectRoot, ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 4D.pdf');
outputFig = fullfile(projectRoot, ...
    'Spectral_Analysis_Report/Figures/Final draft/Figure 4D.fig');

source = load(sourceFile,'output');
trajectory = source.output;
snapshotIndices = [1 5];
snapshotTimesMs = [trajectory.AllSnapshots(snapshotIndices).TimeMs];
assert(isequal(snapshotTimesMs,[3 150]), ...
    'Expected cached snapshots at 3 ms and 150 ms.');

finiteRows = cell(2,1);
for row = 1:2
    cacheFile = fullfile(artifactRoot,'data', ...
        sprintf('figure_04_row_%02d.mat',snapshotIndices(row)));
    loaded = load(cacheFile,'result');
    finiteRows{row} = loaded.result;
    assert(abs(finiteRows{row}.StartTimeMs-snapshotTimesMs(row))<1e-12, ...
        'Finite-time cache does not match the selected snapshot.');
end

firingMaps = cell(2,1);
inputMaps = cell(2,1);
responseMaps = cell(2,1);
for row = 1:2
    firingMaps{row} = trajectory.AllSnapshots(snapshotIndices(row)).FiringMaps{4};
    inputMaps{row} = finiteRows{row}.FiniteTimeOptimalInputEMap;
    responseMaps{row} = finiteRows{row}.FiniteTimeOptimalOutputEMap;
end

firingLimits = shared_limits(firingMaps,false);
inputLimits = shared_limits(inputMaps,true);
responseLimits = shared_limits(responseMaps,true);

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;

fig = figure('Visible','on','Color','w','Units','inches', ...
    'Position',[0.4 0.4 32 14.5],'Renderer','painters');

% HC-norm trajectory panel.
axTrajectory = axes(fig,'Position',[0.035 0.335 0.215 0.365]);
axTrajectory.Toolbar.Visible = 'off';
timeMask = trajectory.Times<=180;
plot(axTrajectory,trajectory.Times(timeMask),trajectory.HCnorm(timeMask), ...
    'Color',neuroColors(6,:),'LineWidth',4.2);
hold(axTrajectory,'on');
markerColors = [neuroColors(1,:); neuroColors(5,:)];
for row = 1:2
    timeMs = snapshotTimesMs(row);
    hcValue = interp1(trajectory.Times,trajectory.HCnorm,timeMs,'linear');
    xline(axTrajectory,timeMs,'--','Color',markerColors(row,:), ...
        'LineWidth',3.0);
    plot(axTrajectory,timeMs,hcValue,'s','MarkerSize',13, ...
        'MarkerFaceColor','w','MarkerEdgeColor',markerColors(row,:), ...
        'LineWidth',3.0);
    if row==1
        text(axTrajectory,timeMs+4,0.91*max(trajectory.HCnorm(timeMask)), ...
            '3 ms','Color',markerColors(row,:),'FontSize',26, ...
            'FontWeight','bold','HorizontalAlignment','left');
    else
        text(axTrajectory,timeMs-4,0.91*max(trajectory.HCnorm(timeMask)), ...
            '150 ms','Color',markerColors(row,:),'FontSize',26, ...
            'FontWeight','bold','HorizontalAlignment','right');
    end
end
xlim(axTrajectory,[0 180]);
ylim(axTrajectory,[0 1.06*max(trajectory.HCnorm(timeMask))]);
xlabel(axTrajectory,'Time (ms)','FontSize',30,'FontWeight','bold');
ylabel(axTrajectory,'HC norm (sp/s)','FontSize',30,'FontWeight','bold');
set(axTrajectory,'FontSize',26,'LineWidth',2.0,'TickDir','out', ...
    'Box','off','XTick',0:30:180,'Layer','top');
grid(axTrajectory,'on');
axTrajectory.GridColor = [0.82 0.82 0.82];
axTrajectory.GridAlpha = 0.55;

columnTitles = {'E firing-rate map', ...
    {'Transient input','direction'}, ...
    {'Transient response','direction'}};
maps = {firingMaps,inputMaps,responseMaps};
limits = {firingLimits,inputLimits,responseLimits};

columnX = [0.310 0.545 0.780];
rowY = [0.575 0.100];
mapWidth = 0.150;
mapHeight = 0.331;
colorbarWidth = 0.010;

for row = 1:2
    for column = 1:3
        ax = axes(fig,'Position',[columnX(column) rowY(row) mapWidth mapHeight]);
        ax.Toolbar.Visible = 'off';
        imagesc(ax,maps{column}{row});
        axis(ax,'image');
        set(ax,'YDir','normal','FontSize',26,'LineWidth',1.8, ...
            'TickDir','out','XTick',[1 10 20 30 40], ...
            'YTick',[1 10 20 30 40]);
        colormap(ax,jet(256));
        clim(ax,limits{column});
        add_hc_grid(ax,size(maps{column}{row}));

        if row==1
            titleSize = 26;
            if column > 1
                titleSize = 22;
            end
            title(ax,columnTitles{column},'FontSize',titleSize, ...
                'FontWeight','bold');
        else
            xlabel(ax,'Map column (pixel)','FontSize',27,'FontWeight','bold');
        end
        if column==1
            ylabel(ax,'Map row (pixel)','FontSize',27,'FontWeight','bold');
            text(ax,-0.34,0.5,sprintf('t = %d ms',snapshotTimesMs(row)), ...
                'Units','normalized','Rotation',90, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'FontSize',29,'FontWeight','bold','Clipping','off');
        end

        if row==1
            cb = colorbar(ax,'eastoutside');
            ax.Position = [columnX(column) rowY(row) mapWidth mapHeight];
            cb.Position = [columnX(column)+mapWidth+0.014, ...
                rowY(row),colorbarWidth,mapHeight];
            cb.FontSize = 26;
            cb.LineWidth = 1.5;
            cb.Ruler.Exponent = 0;
            if max(abs(limits{column}))<0.01
                cb.Ruler.TickLabelFormat = '%.3f';
            end
            cb.Title.String = 'sp/s';
            cb.Title.FontSize = 27;
            cb.Title.FontWeight = 'bold';
        end
    end
end

savefig(fig,outputFig);
exportgraphics(fig,outputPdf,'ContentType','image','Resolution',600);
close(fig);
movefile(outputPdf,finalPdf,'f');
fprintf('Saved %s\n',finalPdf);
end

function limits = shared_limits(maps,isSigned)
values = vertcat(maps{1}(:),maps{2}(:));
if isSigned
    maximum = max(abs(values));
    limits = maximum*[-1 1];
else
    limits = [min(values) max(values)];
end
if diff(limits)<=eps
    limits = mean(limits)+[-1 1]*1e-6;
end
end

function add_hc_grid(ax,mapSize)
hold(ax,'on');
for value = 10.5:10:(mapSize(2)-0.5)
    xline(ax,value,'-','Color',[0.48 0.48 0.48],'LineWidth',0.75);
end
for value = 10.5:10:(mapSize(1)-0.5)
    yline(ax,value,'-','Color',[0.48 0.48 0.48],'LineWidth',0.75);
end
hold(ax,'off');
end
