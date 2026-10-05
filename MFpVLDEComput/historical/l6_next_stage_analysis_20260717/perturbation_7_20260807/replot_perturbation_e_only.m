function replot_perturbation_e_only(dataFile,outputPdf,figureNumberOverride)
% Replot one saved perturbation experiment using only the E population.

arguments
    dataFile (1,:) char
    outputPdf (1,:) char
    figureNumberOverride (1,:) char = ''
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
times = result.Times;
hc = result.HCnorm;
singularAlignment = result.SingularAlignment;
returnAlignment = result.ReturnAlignment;

firingLimits = local_pair_limits(snapshots,'FiringMaps',false);
returnMaximum = local_pair_limits(snapshots,'ReturnDirectionMaps',true);
singularLimits = local_pair_limits(snapshots,'SingularInputMaps',false);

snapshotCount = numel(snapshots);
fig = figure('Visible','off','Color','w', ...
    'Position',[20 20 1800 650+300*snapshotCount]);
layout = tiledlayout(fig,1+snapshotCount,3, ...
    'TileSpacing','compact','Padding','loose');
layout.OuterPosition = [0 0 1 0.94];

ax = nexttile(layout,1);
plot(ax,times,hc,'k-','LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'HC norm from fixed point'); grid(ax,'on');
title(ax,'Nonlinear ODE perturbation amplitude','FontWeight','bold');

ax = nexttile(layout,2);
plot(ax,times,singularAlignment,'Color',[0.85 0.33 0.10],'LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'Alignment'); ylim(ax,[0 1]); grid(ax,'on');
title(ax,'Alignment with fixed-point transient output','FontWeight','bold');

ax = nexttile(layout,3);
plot(ax,times,returnAlignment,'Color',[0 0.45 0.74],'LineWidth',1.8); hold(ax,'on');
local_time_lines(ax,snapshots,summary);
xlabel(ax,'Time (ms)'); ylabel(ax,'Cosine alignment'); ylim(ax,[-1 1]); grid(ax,'on');
title(ax,'ODE velocity alignment with fixed-point return direction','FontWeight','bold');

columnTitles = {'E firing-rate map (Hz)', ...
    'E direction toward stable fixed point','E top right singular vector'};
for snapshotIndex = 1:snapshotCount
    maps = {snapshots(snapshotIndex).FiringMaps{4}, ...
        snapshots(snapshotIndex).ReturnDirectionMaps{4}, ...
        snapshots(snapshotIndex).SingularInputMaps{4}};
    limits = {firingLimits, returnMaximum*[-1 1], ...
        singularLimits};
    for columnIndex = 1:3
        if local_is_fixed_point_snapshot(snapshots(snapshotIndex)) && columnIndex>1
            limits{columnIndex} = local_map_limits(maps{columnIndex},columnIndex==2);
        end
        ax = nexttile(layout,3+3*(snapshotIndex-1)+columnIndex);
        imagesc(ax,maps{columnIndex});
        axis(ax,'image');
        set(ax,'XTick',[],'YTick',[],'YDir','normal');
        clim(ax,limits{columnIndex});
        local_hc_grid(ax,size(maps{columnIndex}));
        colorbar(ax,'eastoutside');
        if snapshotIndex==1
            title(ax,columnTitles{columnIndex},'FontWeight','bold','FontSize',12);
        end
        if columnIndex==1
            ylabel(ax,sprintf('%s, E, t=%.1f ms', ...
                local_phase_name(snapshots(snapshotIndex),summary), ...
                snapshots(snapshotIndex).TimeMs), ...
                'Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold','FontSize',11);
        end
    end
end

colormap(fig,jet(256));
figureTitle = sprintf([ ...
    '%s: Contrast 100, orientation 0.0 deg - %s; E population only\n' ...
    'Initial HC=%.3f, w_6=%+.2f, L6 gain=%.2fx, max Re lambda=%.6f, ' ...
    'sigma_1=%.4f, tau=%.3f ms'], ...
    figureNumber,caseLabel,result.InitialHCnorm,result.L6Weight,1-result.L6Weight, ...
    real(result.FixedPointLeadingEigenvalue),result.FixedPointTopSingularValue, ...
    result.TauMs);
titleAxis = axes(fig,'Position',[0.03 0.925 0.94 0.065]);
axis(titleAxis,'off');
text(titleAxis,0.5,0.5,figureTitle,'Units','normalized', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontWeight','bold','FontSize',12,'Interpreter','tex');
set(findall(fig,'Type','axes'),'FontSize',10,'LineWidth',0.8);
exportgraphics(fig,outputPdf,'ContentType','image','Resolution',300);
close(fig);
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
        color = [0.85 0.33 0.10];
        style = '--';
    elseif abs(timeMs-summary.convergenceTimeMs(1))<1e-9
        color = [0 0.45 0.74];
        style = '--';
    else
        color = [0.93 0.23 0.23];
        style = '-';
    end
    xline(axisHandle,timeMs,style,'LineWidth',1.3,'Color',color);
end
end

function name = local_phase_name(snapshot,summary)
timeMs = snapshot.TimeMs;
if isfield(snapshot,'PhaseName') && strlength(string(snapshot.PhaseName))>0
    name = char(snapshot.PhaseName);
    return;
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
