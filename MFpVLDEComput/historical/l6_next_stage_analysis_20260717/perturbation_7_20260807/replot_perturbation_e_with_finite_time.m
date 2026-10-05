function replot_perturbation_e_with_finite_time(dataFile,finiteTimeRoot,outputPdf,figureNumberOverride,baselineSingularOutputEMap,denseCurveRoot)
% Add trajectory-aware finite-time transient response maps to one figure.

arguments
    dataFile (1,:) char
    finiteTimeRoot (1,:) char
    outputPdf (1,:) char
    figureNumberOverride (1,:) char = ''
    baselineSingularOutputEMap double = []
    denseCurveRoot (1,:) char = ''
end

loaded = load(dataFile,'output');
result = loaded.output;
if isfield(result,'AllSnapshots')
    snapshots = result.AllSnapshots;
else
    snapshots = result.Snapshots;
end
summary = result.Summary;
figureNumber = char(summary.figureNumber(1));
if ~isempty(figureNumberOverride)
    figureNumber = figureNumberOverride;
end
caseLabel = char(summary.caseLabel(1));
figureIndex = local_figure_index(figureNumber);

finiteRows = cell(numel(snapshots),1);
for row = 1:numel(snapshots)
    loadedRow = load(fullfile(finiteTimeRoot, ...
        sprintf('figure_%02d_row_%02d.mat',figureIndex,row)),'result');
    finiteRows{row} = loadedRow.result;
end

times = result.Times;
hc = result.HCnorm;
singularAlignment = result.SingularAlignment;
returnAlignment = result.ReturnAlignment;
denseLoaded = load(fullfile(denseCurveRoot, ...
    sprintf('figure_%02d_curve.mat',figureIndex)),'curve');
denseCurve = denseLoaded.curve;
firingLimits = local_pair_limits(snapshots,'FiringMaps',false);
returnMaximum = local_pair_limits(snapshots,'ReturnDirectionMaps',true);

snapshotCount = numel(snapshots);
fig = figure('Visible','off','Color','w', ...
    'Position',[20 20 2250 650+300*snapshotCount]);
layout = tiledlayout(fig,1+snapshotCount,12, ...
    'TileSpacing','compact','Padding','loose');
layout.OuterPosition = [0 0 1 0.94];

ax = nexttile(layout,1,[1 4]);
plot(ax,times,hc,'k-','LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'HC norm from fixed point'); grid(ax,'on');
title(ax,'Nonlinear ODE perturbation amplitude','FontWeight','bold');

ax = nexttile(layout,5,[1 4]);
plot(ax,times,singularAlignment,'Color',[0.85 0.33 0.10],'LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'Alignment'); ylim(ax,[0 1]); grid(ax,'on');
title(ax,'Alignment with fixed-point transient output','FontWeight','bold');

ax = nexttile(layout,9,[1 4]);
plot(ax,denseCurve.timesMs,denseCurve.alignment,'-', ...
    'Color',[0 0.45 0.74],'LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'Absolute cosine'); ylim(ax,[0 1]);
xlim(ax,[0 500]); grid(ax,'on');
title(ax,['Finite-time response alignment with baseline top left ' ...
    'singular vector'],'FontWeight','bold');

columnTitles = {'E firing-rate map (Hz)', ...
    'E direction toward stable fixed point', ...
    'E instantaneous top right singular vector', ...
    'E finite-time transient response direction'};
for row = 1:snapshotCount
    maps = {snapshots(row).FiringMaps{4}, ...
        snapshots(row).ReturnDirectionMaps{4}, ...
        snapshots(row).SingularInputMaps{4}, ...
        finiteRows{row}.FiniteTimeOptimalOutputEMap};
    limits = {firingLimits,returnMaximum*[-1 1], ...
        local_map_limits(maps{3},true),local_map_limits(maps{4},true)};
    for column = 1:4
        if local_is_fixed_point_snapshot(snapshots(row)) && column==2
            limits{column} = local_map_limits(maps{column},true);
        end
        tileIndex = 12*row+1+3*(column-1);
        ax = nexttile(layout,tileIndex,[1 3]);
        imagesc(ax,maps{column});
        axis(ax,'image');
        set(ax,'XTick',[],'YTick',[],'YDir','normal');
        clim(ax,limits{column});
        local_hc_grid(ax,size(maps{column}));
        colorbar(ax,'eastoutside');
        if row==1
            title(ax,columnTitles{column},'FontWeight','bold','FontSize',11);
        end
        if column==1
            ylabel(ax,sprintf('%s, E, t=%.1f ms', ...
                local_phase_name(snapshots(row),summary),snapshots(row).TimeMs), ...
                'Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold','FontSize',10);
        elseif column==4
            xlabel(ax,sprintf('Tpeak=%.1f ms; gain=%.2f', ...
                finiteRows{row}.PeakHorizonMs,finiteRows{row}.PeakGain), ...
                'FontSize',8,'Interpreter','none');
        end
    end
end

colormap(fig,jet(256));
figureTitle = sprintf([ ...
    '%s: Contrast 100, orientation 0.0 deg - %s; E population only\n' ...
    'Initial HC=%.3f, w6=%+.2f, L6 gain=%.2fx, max Re lambda=%.6f, ' ...
    'sigma1=%.4f, tau=%.3f ms'], ...
    figureNumber,caseLabel,result.InitialHCnorm,result.L6Weight,1-result.L6Weight, ...
    real(result.FixedPointLeadingEigenvalue),result.FixedPointTopSingularValue, ...
    result.TauMs);
titleAxis = axes(fig,'Position',[0.03 0.925 0.94 0.065]);
axis(titleAxis,'off');
text(titleAxis,0.5,0.5,figureTitle,'Units','normalized', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',12,'Interpreter','none');
set(findall(fig,'Type','axes'),'FontSize',9,'LineWidth',0.8);
exportgraphics(fig,outputPdf,'ContentType','image','Resolution',300);
close(fig);
end

function index = local_figure_index(figureNumber)
numbers = {'7.1','7.2','7.3','7.4','7.5','7.6', ...
    '7.7','7.8','7.9','7.10','7.11','7.12'};
index = find(strcmp(numbers,figureNumber),1);
if isempty(index)
    error('Unknown figure number %s.',figureNumber);
end
end

function tf = local_is_fixed_point_snapshot(snapshot)
tf = isfield(snapshot,'PhaseName') && ...
    strcmp(string(snapshot.PhaseName),"Fixed point");
end

function limits = local_map_limits(map,isSigned)
if isSigned
    maximum = max(abs(map(:)));
    limits = maximum*[-1 1];
else
    limits = [min(map(:)) max(map(:))];
end
if diff(limits)<=eps
    center = mean(limits);
    limits = center+[-1 1]*1e-6;
end
end

function limits = local_pair_limits(snapshots,fieldName,isSigned)
values = [];
for index = 1:numel(snapshots)
    values = [values; snapshots(index).(fieldName){4}(:)]; %#ok<AGROW>
end
if isSigned
    limits = max(max(abs(values)),eps);
else
    limits = [min(values) max(values)];
    if diff(limits)<=0
        limits = limits+[-1 1]*1e-6;
    end
end
end

function local_time_lines(axisHandle,snapshots,summary)
for index = 1:numel(snapshots)
    timeMs = snapshots(index).TimeMs;
    if abs(timeMs-summary.transientTimeMs(1))<1e-9
        color = [0.85 0.33 0.10]; style = '--';
    elseif abs(timeMs-summary.convergenceTimeMs(1))<1e-9
        color = [0 0.45 0.74]; style = '--';
    else
        color = [0.93 0.23 0.23]; style = '-';
    end
    xline(axisHandle,timeMs,style,'LineWidth',1.3,'Color',color);
end
end

function name = local_phase_name(snapshot,summary)
timeMs = snapshot.TimeMs;
if isfield(snapshot,'PhaseName') && strlength(string(snapshot.PhaseName))>0
    name = char(snapshot.PhaseName); return;
end
if abs(timeMs-summary.transientTimeMs(1))<1e-9
    name = 'Ongoing transient';
elseif abs(timeMs-summary.convergenceTimeMs(1))<1e-9
    if summary.convergenceReturnAlignment(1)>0 && ...
            summary.convergenceHC(1)<summary.actualInitialHC(1)
        name = 'Early return toward fixed point';
    else
        name = 'Late non-returning state';
    end
else
    name = 'Additional snapshot';
end
end

function local_hc_grid(axisHandle,mapSize)
hold(axisHandle,'on');
for value = 10.5:10:(mapSize(2)-0.5)
    xline(axisHandle,value,'-','Color',[0.55 0.55 0.55],'LineWidth',0.35);
end
for value = 10.5:10:(mapSize(1)-0.5)
    yline(axisHandle,value,'-','Color',[0.55 0.55 0.55],'LineWidth',0.35);
end
hold(axisHandle,'off');
end
