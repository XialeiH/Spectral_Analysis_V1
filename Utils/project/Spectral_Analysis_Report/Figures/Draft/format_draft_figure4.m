function format_draft_figure4(fig)
% Shared typography and physical spacing for the six requested Draft figures.
axesList = findall(fig,'Type','axes');
positions = zeros(numel(axesList),4);
for k=1:numel(axesList)
    axesList(k).Units='inches';
    positions(k,:)=axesList(k).Position;
    axesList(k).XTick=axesList(k).XTick;
    for ruler=axesList(k).YAxis(:).', ruler.TickValues=ruler.TickValues; end
end
set(findall(fig,'Type','text'),'FontName','Arial','FontSize',51,'FontWeight','bold');
for ax = axesList(:).'
    ax.Units = 'inches';
    set(ax,'FontName','Arial','FontSize',41,'FontWeight','normal');
    if ~isempty(findall(ax,'Type','image'))
        set(ax,'XTick',[],'YTick',[]);
    elseif ~isempty(ax.XTickLabel) && any(~cellfun('isempty',regexp(cellstr(ax.XTickLabel),'[A-Za-z]','once')))
        labels = cellstr(ax.XTickLabel);
        labels = cellfun(@(s) ['{\bf ' strrep(s,'6','{\rm6}') '}'],labels,'UniformOutput',false);
        ax.XTickLabel = labels;
        ax.TickLabelInterpreter = 'tex';
    end
    set([ax.Title ax.XLabel],'FontSize',51,'FontWeight','bold');
    if isequal(ax.Title.String,'Right E Eigencluster')
        ax.Title.String={'Right E';'Eigencluster'};
    elseif isequal(ax.Title.String,'Left E Eigencluster')
        ax.Title.String={'Left E';'Eigencluster'};
    elseif isequal(ax.Title.String,'Right Singular Mode')
        ax.Title.String={'Right Singular';'Mode'};
    elseif isequal(ax.Title.String,'Left Singular Mode')
        ax.Title.String={'Left Singular';'Mode'};
    end
    for ruler = ax.YAxis(:).'
        set(ruler.Label,'FontSize',51,'FontWeight','bold');
    end
    if isequal(ax.YLabel.String,'Eigencluster similarity')
        ax.YLabel.String = {'Eigencluster';'similarity'};
    elseif isequal(ax.YLabel.String,'Normalized mode envelope')
        ax.YLabel.String = {'Normalized mode';'envelope'};
    end
end
for label=findall(fig,'Type','text').'
    if ischar(label.String)
        label.String=regexprep(label.String,'t = (\d+) ms','t = {\\rm$1} ms');
        label.String=strrep(label.String,'J_6','J_{\rm6}');
    end
end

% Enlarge legend samples independently of the scientific curves.
legends = findall(fig,'Type','legend');
for lg = legends(:).'
    ax = lg.Axes;
    entries = lg.PlotChildren;
    labels = cellstr(lg.String);
    proxies = gobjects(size(entries));
    hold(ax,'on');
    for k = 1:numel(entries)
        h = entries(k);
        proxies(k) = plot(ax,nan,nan,'LineStyle',h.LineStyle,'Color',h.Color, ...
            'LineWidth',5,'Marker',h.Marker,'MarkerSize',h.MarkerSize, ...
            'MarkerFaceColor',h.MarkerFaceColor,'MarkerEdgeColor',h.MarkerEdgeColor, ...
            'Tag','LegendOnly');
        labels{k} = ['  ' strrep(labels{k},'L6','L{\rm6}')];
    end
    hold(ax,'off');
    lg = legend(ax,proxies,labels,'FontName','Arial','FontSize',41, ...
        'FontWeight','bold','Interpreter','tex','Box','on','Color','white', ...
        'AutoUpdate','off','NumColumns',1,'Location','northeastoutside');
    lg.ItemTokenSize = [64 18];
    lg.Units = 'inches';
    drawnow;
    p = lg.Position;
    % Leave room for the right-side axis title on the dual-axis trajectory.
    extra = 0.5 + 2.0*(numel(ax.YAxis)>1);
    p(1:2) = [ax.Position(1)+ax.Position(3)+extra, ...
        ax.Position(2)+ax.Position(4)-p(4)];
    lg.Position = p;
end
for k=1:numel(axesList), axesList(k).Position=positions(k,:); end

% Measure the 0.4-inch label gap from the tick-text extent, not the axis line.
drawnow;
for ax = axesList(:).'
    labels = [ax.XLabel; arrayfun(@(r) r.Label,ax.YAxis(:))];
    strings = get(labels,'String');
    for k = 1:numel(labels), labels(k).String = ''; end
    drawnow;
    inset = ax.TightInset;
    for k = 1:numel(labels), labels(k).String = strings{k}; end
    set(labels,'Units','inches');
    drawnow;
    for k = 1:numel(labels)
        label = labels(k);
        if isempty(label.String), continue; end
        e = label.Extent; p = label.Position;
        if k == 1
            p(2) = p(2)-inset(2)-0.4-e(2)-e(4);
        elseif (k==3) || (numel(ax.YAxis)==1 && strcmp(ax.YAxisLocation,'right'))
            p(1) = p(1)+ax.Position(3)+inset(3)+0.4-e(1);
        else
            p(1) = p(1)-inset(1)-0.4-e(1)-e(3);
        end
        label.Position = p;
    end
end

bars = findall(fig,'Type','colorbar');
imageAxes = axesList(arrayfun(@(a) ~isempty(findall(a,'Type','image')),axesList));
for cb = bars(:).'
    cb.Units = 'inches';
    old = cb.Position;
    distances = arrayfun(@(a) norm(old(1:2)- ...
        [a.Position(1)+a.Position(3),a.Position(2)]),imageAxes);
    [~,index] = min(distances);
    ax = imageAxes(index);
    cb.Position = [ax.Position(1)+ax.Position(3)+0.3 ax.Position(2) 30/72 ax.Position(4)];
    set(cb,'FontName','Arial','FontSize',41,'FontWeight','normal');
end

% Retain panel sizes; grow the canvas only to accommodate larger text/legends.
drawnow;
canvas = fig.Position(3:4);
for lg = findall(fig,'Type','legend').'
    canvas = max(canvas,lg.Position(1:2)+lg.Position(3:4)+[0.3 0.3]);
end
fig.Position(3:4) = canvas;
for k=1:numel(axesList), axesList(k).Position=positions(k,:); end
drawnow;
for lg = findall(fig,'Type','legend').'
    ax=lg.Axes;
    p=lg.Position;
    p(1:2)=[ax.Position(1)+ax.Position(3)+0.5+2.0*(numel(ax.YAxis)>1), ...
        ax.Position(2)+ax.Position(4)-p(4)];
    lg.Position=p;
    canvas=max(canvas,p(1:2)+p(3:4)+0.3);
end
fig.Position(3:4)=canvas;
for k=1:numel(axesList), axesList(k).Position=positions(k,:); end
if ~isempty(bars)
    textAxes = axes(fig,'Units','inches','Position',[0 0 canvas], ...
        'XLim',[0 canvas(1)],'YLim',[0 canvas(2)],'Visible','off', ...
        'Color','none','Tag','ColorbarText');
    for cb = bars(:).'
        p = cb.Position; ticks = cb.Ticks;
        cb.TickLabels = repmat({''},size(ticks));
        cb.Title.Visible = 'off';
        for tick = ticks
            y = p(2)+(tick-cb.Limits(1))/diff(cb.Limits)*p(4);
            text(textAxes,p(1)+p(3)+0.06,y,sprintf('%g',tick), ...
                'FontName','Arial','FontSize',41,'FontWeight','normal', ...
                'HorizontalAlignment','left','VerticalAlignment','middle','Clipping','off');
        end
        text(textAxes,p(1)+p(3)/2,p(2)+p(4)+0.1,cb.Title.String, ...
            'FontName','Arial','FontSize',51,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','off');
    end
end
drawnow;
% Crop the page physically, without exportgraphics scaling text to old bounds.
for lg=findall(fig,'Type','legend').'
    ax=lg.Axes;
    if isequal(ax.YLabel.String,{'Normalized mode';'envelope'})
        p=lg.Position;
        p(1:2)=[ax.Position(1)+(ax.Position(3)-p(3))/2, ...
            ax.Position(2)+ax.Position(4)+0.4];
        lg.Position=p;
    end
end
drawnow;
allAxes=findall(fig,'Type','axes');
low=[inf inf]; high=[-inf -inf];
for ax=allAxes(:).'
    ax.Units='inches'; p=ax.Position;
    if strcmp(ax.Visible,'on')
        inset=ax.TightInset;
        low=min(low,p(1:2)-inset(1:2));
        high=max(high,p(1:2)+p(3:4)+inset(3:4));
    end
    for t=findall(ax,'Type','text').'
        if strcmp(t.Visible,'off') || isempty(t.String), continue; end
        t.Units='inches'; e=t.Extent;
        low=min(low,p(1:2)+e(1:2));
        high=max(high,p(1:2)+e(1:2)+e(3:4));
    end
end
allLegends=findall(fig,'Type','legend');
for lg=allLegends(:).'
    lg.Units='inches'; p=lg.Position;
    low=min(low,p(1:2));high=max(high,p(1:2)+p(3:4));
end
shift=[0.25 0.25]-low;
canvas=high-low+0.5;
axPositions=arrayfun(@(a) a.Position,allAxes,'UniformOutput',false);
barPositions=arrayfun(@(a) a.Position,bars,'UniformOutput',false);
legendPositions=arrayfun(@(a) a.Position,allLegends,'UniformOutput',false);
fig.Position=[0.2 0.2 canvas];
fig.Position=[0.2 0.2 canvas];
for k=1:numel(allAxes), allAxes(k).Position=axPositions{k}+[shift 0 0]; end
for k=1:numel(bars), bars(k).Position=barPositions{k}+[shift 0 0]; end
for k=1:numel(allLegends), allLegends(k).Position=legendPositions{k}+[shift 0 0]; end
fig.PaperUnits='inches';fig.PaperSize=canvas;fig.PaperPosition=[0 0 canvas];
fig.PaperPositionMode='manual';
drawnow;
fprintf('Final page %.4f x %.4f inches; actual figure %.4f x %.4f.\n', ...
    canvas,fig.Position(3:4));
end
