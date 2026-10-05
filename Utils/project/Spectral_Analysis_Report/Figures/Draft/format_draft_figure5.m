function format_draft_figure5(fig,fontIncrement)
% Format the requested Figure 5 panels without modifying their plotted data.
if nargin<2,fontIncrement=0;end
labelFont=51+fontIncrement;
numberFont=41+fontIncrement;
axs=findall(fig,'Type','axes');
positions=zeros(numel(axs),4);
for k=1:numel(axs)
    ax=axs(k);ax.Units='inches';positions(k,:)=ax.Position;
    ax.XTick=ax.XTick;ax.YTick=ax.YTick;
    set(ax,'FontName','Arial','FontSize',numberFont,'FontWeight','normal');
    if ~isempty(findall(ax,'Type','image')) && ~strcmp(ax.Tag,'PreserveAxisNumbers')
        ax.XTick=[];ax.YTick=[];
    end
    for t=findall(ax,'Type','text').'
        set(t,'FontName','Arial','FontSize',labelFont,'FontWeight','bold');
        if ischar(t.String)
            if ~isempty(regexp(t.String,'^\d+(\.\d+)?$','once'))
                set(t,'FontSize',numberFont,'FontWeight','normal');
            else
                t.String=regexprep(t.String,'(Condition\s+)([123])','$1{\\rm$2}');
            end
        end
    end
end

for lg=findall(fig,'Type','legend').'
    ax=lg.Axes;entries=lg.PlotChildren;labels=cellstr(lg.String);
    proxies=gobjects(size(entries));hold(ax,'on');
    for k=1:numel(entries)
        h=entries(k);style='-';color=[0 0 0];marker='none';face='none';edge='auto';ms=12;
        if isprop(h,'LineStyle'),style=h.LineStyle;end
        if isprop(h,'Color'),color=h.Color;
        elseif isprop(h,'LineColor'),color=h.LineColor;
        elseif isprop(h,'CData') && numel(h.CData)==3,color=h.CData;
        end
        if isprop(h,'Marker'),marker=h.Marker;end
        if isprop(h,'MarkerSize'),ms=max(12,h.MarkerSize);end
        if isprop(h,'MarkerFaceColor'),face=h.MarkerFaceColor;end
        if isprop(h,'MarkerEdgeColor'),edge=h.MarkerEdgeColor;end
        if strcmp(h.Type,'scatter')
            style='none';face=color;edge=color;ms=max(14,sqrt(h.SizeData(1)));
        end
        if strcmp(h.Type,'patch')
            style='none';marker='s';face=h.FaceColor;edge=face;
        end
        if ischar(color),color=[0.4 0.4 0.4];end
        proxies(k)=plot(ax,nan,nan,'LineStyle',style,'Color',color, ...
            'LineWidth',5,'Marker',marker,'MarkerSize',ms, ...
            'MarkerFaceColor',face,'MarkerEdgeColor',edge,'Tag','LegendOnly');
        labels{k}=['  ' regexprep(strtrim(labels{k}),'(\d+(?:\.\d+)?)','{\\rm$1}')];
        labels{k}=strrep(labels{k},'Orientation Tuning Curve Error (',sprintf('Orientation Tuning Curve Error\n('));
    end
    hold(ax,'off');
    ncols=1;
    if numel(labels)==5
        ncols=2;
        if strcmp(ax.Tag,'PreserveAxisNumbers'),ncols=3;end
    end
    lg=legend(ax,proxies,labels,'FontName','Arial','FontSize',numberFont, ...
        'FontWeight','bold','Interpreter','tex','Box','on','Color','white', ...
        'AutoUpdate','off','NumColumns',ncols,'Location','northoutside');
    lg.ItemTokenSize=[64 18];lg.Units='inches';
end
for k=1:numel(axs),axs(k).Position=positions(k,:);end
drawnow;
for ax=axs(:).'
    labels=[ax.XLabel ax.YLabel];strings=get(labels,'String');
    set(labels,'String','');drawnow;inset=ax.TightInset;
    for k=1:2,labels(k).String=strings{k};end
    set(labels,'Units','inches');drawnow;
    for k=1:2
        t=labels(k);if isempty(t.String),continue;end
        e=t.Extent;p=t.Position;
        if k==1,p(2)=p(2)-inset(2)-0.4-e(2)-e(4);
        else,p(1)=p(1)-inset(1)-0.4-e(1)-e(3);end
        t.Position=p;
    end
end

bars=findall(fig,'Type','colorbar');
imageAxes=axs(arrayfun(@(a) ~isempty(findall(a,'Type','image')),axs));
canvas=fig.Position(3:4);
textAxes=axes(fig,'Units','inches','Position',[0 0 canvas], ...
    'XLim',[0 canvas(1)],'YLim',[0 canvas(2)],'Visible','off', ...
    'Color','none','Tag','ColorbarText');
for cb=bars(:).'
    cb.Units='inches';p=cb.Position;
    d=arrayfun(@(a) norm(p(1:2)-[a.Position(1)+a.Position(3) a.Position(2)]),imageAxes);
    [~,j]=min(d);ax=imageAxes(j);
    cb.Position=[ax.Position(1)+ax.Position(3)+0.3 ax.Position(2) 30/72 ax.Position(4)];
    cb.AxisLocation='out';cb.FontSize=numberFont;cb.FontWeight='normal';
    p=cb.Position;ticks=cb.Ticks;cb.TickLabels=repmat({''},size(ticks));
    cb.Title.Visible='off';cb.Label.Visible='off';
    right=p(1)+p(3);
    for tick=ticks
        if strcmp(ax.ColorScale,'log')
            fraction=log(tick/cb.Limits(1))/log(cb.Limits(2)/cb.Limits(1));
            label=sprintf('10^{%g}',log10(tick));
        else
            fraction=(tick-cb.Limits(1))/diff(cb.Limits);label=sprintf('%g',tick);
        end
        t=text(textAxes,p(1)+p(3)+0.06,p(2)+fraction*p(4),label, ...
            'FontName','Arial','FontSize',numberFont,'FontWeight','normal', ...
            'HorizontalAlignment','left','VerticalAlignment','middle','Clipping','off');
        t.Units='inches';drawnow;right=max(right,t.Extent(1)+t.Extent(3));
    end
    text(textAxes,p(1)+p(3)/2,p(2)+p(4)+0.1,cb.Title.String, ...
        'FontName','Arial','FontSize',labelFont,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','off');
    if ~isempty(cb.Label.String)
        text(textAxes,right+0.4,p(2)+p(4)/2,cb.Label.String, ...
            'Rotation',90,'FontName','Arial','FontSize',labelFont,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment','top','Clipping','off');
    end
end
drawnow;
for lg=findall(fig,'Type','legend').'
    p=lg.Position;
    left=min(positions(:,1));right=max(positions(:,1)+positions(:,3));
    top=max(positions(:,2)+positions(:,4));
    lg.Position=[(left+right-p(3))/2 top+0.45 p(3:4)];
end
drawnow;
% Fit physical page bounds so MATLAB does not rescale fonts during export.
allAxes=findall(fig,'Type','axes');low=[inf inf];high=[-inf -inf];
for ax=allAxes(:).'
    ax.Units='inches';p=ax.Position;
    if strcmp(ax.Visible,'on')
        inset=ax.TightInset;low=min(low,p(1:2)-inset(1:2));
        high=max(high,p(1:2)+p(3:4)+inset(3:4));
    end
    for t=findall(ax,'Type','text').'
        if strcmp(t.Visible,'off') || isempty(t.String),continue;end
        t.Units='inches';e=t.Extent;
        low=min(low,p(1:2)+e(1:2));high=max(high,p(1:2)+e(1:2)+e(3:4));
    end
end
legends=findall(fig,'Type','legend');
for lg=legends(:).'
    p=lg.Position;low=min(low,p(1:2));high=max(high,p(1:2)+p(3:4));
end
shift=0.25-low;canvas=high-low+0.5;
ap=arrayfun(@(a) a.Position,allAxes,'UniformOutput',false);
bp=arrayfun(@(a) a.Position,bars,'UniformOutput',false);
lp=arrayfun(@(a) a.Position,legends,'UniformOutput',false);
fig.Position=[0.2 0.2 canvas];fig.Position=[0.2 0.2 canvas];
for k=1:numel(allAxes),allAxes(k).Position=ap{k}+[shift 0 0];end
for k=1:numel(bars),bars(k).Position=bp{k}+[shift 0 0];end
for k=1:numel(legends),legends(k).Position=lp{k}+[shift 0 0];end
fig.PaperUnits='inches';fig.PaperSize=canvas;fig.PaperPosition=[0 0 canvas];
fig.PaperPositionMode='manual';drawnow;
end
