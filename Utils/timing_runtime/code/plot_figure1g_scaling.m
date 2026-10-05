function plot_figure1g_scaling(trialFile, scalingFile, outputPdf)
% Figure 1G: 4x4 per-iteration timing and square-field end-to-end scaling.
trials = readtable(trialFile,'FileType','text','Delimiter','\t');
scaling = readtable(scalingFile,'FileType','text','Delimiter','\t');
assert(height(trials)==100,'Left panel requires the 100 paired 4x4 trials.');
fieldSizes = [4 6 8 10 20 30 40];
assert(isequal(unique(scaling.FieldHC).',fieldSizes), ...
    'Scaling results must contain all seven requested field sizes.');

colors = [178 24 43;227 74 51;254 227 145;171 217 233;67 147 195;33 102 172]/255;
cgColor = colors(1,:);
dnnColor = colors(6,:);
fig = figure('Color','w','Units','inches','Position',[0.4 0.4 19.0 8.8]);
ax1 = axes(fig,'Position',[0.075 0.22 0.365 0.65]);
draw_four_by_four(ax1,trials.CGSecondsPerIteration,trials.DNNSecondsPerIteration,cgColor,dnnColor);

ax2 = axes(fig,'Position',[0.595 0.22 0.355 0.65]);
[hCG,hDNN] = draw_scaling(ax2,scaling,fieldSizes,cgColor,dnnColor);
lgd = legend(ax2,[hCG hDNN],{'CG','DNN surrogate'}, ...
    'Orientation','horizontal','Box','off','FontSize',26,'Units','normalized');
lgd.Position = [0.70 0.025 0.23 0.055];

exportgraphics(fig,outputPdf,'ContentType','vector');
end

function draw_four_by_four(ax,cgValues,dnnValues,cgColor,dnnColor)
hold(ax,'on');
plot(ax,ones(size(cgValues)),cgValues,'o','LineStyle','none', ...
    'MarkerSize',4.5,'MarkerFaceColor',cgColor,'MarkerEdgeColor','none');
plot(ax,2*ones(size(dnnValues)),dnnValues,'o','LineStyle','none', ...
    'MarkerSize',4.5,'MarkerFaceColor',dnnColor,'MarkerEdgeColor','none');
[cgMean,cgSem] = mean_sem(cgValues);
[dnnMean,dnnSem] = mean_sem(dnnValues);
errorbar(ax,1,cgMean,cgSem,'o','Color',cgColor,'MarkerFaceColor',cgColor, ...
    'MarkerSize',10,'LineWidth',2.8,'CapSize',18);
errorbar(ax,2,dnnMean,dnnSem,'o','Color',dnnColor,'MarkerFaceColor',dnnColor, ...
    'MarkerSize',10,'LineWidth',2.8,'CapSize',18);
set(ax,'YScale','log','XLim',[0.55 2.45],'XTick',[1 2], ...
    'XTickLabel',{'CG','DNN surrogate'},'FontName','Arial','FontSize',26, ...
    'LineWidth',1.8,'TickDir','out','Box','off');
ylabel(ax,'Average time per iteration (s)','FontSize',29);
title(ax,sprintf('4 x 4 HC speedup: %.1fx',median(cgValues./dnnValues)), ...
    'FontSize',28,'FontWeight','normal');
grid(ax,'on'); ax.GridAlpha=0.14;
end

function [hCG,hDNN] = draw_scaling(ax,T,fieldSizes,cgColor,dnnColor)
hold(ax,'on');
x = 1:numel(fieldSizes);
offset = 0.10;
cgMean = zeros(size(x)); cgSem = zeros(size(x));
dnnMean = zeros(size(x)); dnnSem = zeros(size(x));
for k = 1:numel(fieldSizes)
    rows = T.FieldHC==fieldSizes(k);
    cg = T.CGSeconds(rows);
    dnn = T.DNNSeconds(rows);
    plot(ax,(x(k)-offset)*ones(size(cg)),cg,'o','LineStyle','none', ...
        'MarkerSize',5.2,'MarkerFaceColor',cgColor,'MarkerEdgeColor','none');
    plot(ax,(x(k)+offset)*ones(size(dnn)),dnn,'o','LineStyle','none', ...
        'MarkerSize',5.2,'MarkerFaceColor',dnnColor,'MarkerEdgeColor','none');
    [cgMean(k),cgSem(k)] = mean_sem(cg);
    [dnnMean(k),dnnSem(k)] = mean_sem(dnn);
end
hCG = plot(ax,x-offset,cgMean,'-o','Color',cgColor,'MarkerFaceColor',cgColor, ...
    'MarkerSize',8,'LineWidth',2.8);
hDNN = plot(ax,x+offset,dnnMean,'-o','Color',dnnColor,'MarkerFaceColor',dnnColor, ...
    'MarkerSize',8,'LineWidth',2.8);
errorbar(ax,x-offset,cgMean,cgSem,'LineStyle','none','Color',cgColor, ...
    'LineWidth',2.3,'CapSize',12,'HandleVisibility','off');
errorbar(ax,x+offset,dnnMean,dnnSem,'LineStyle','none','Color',dnnColor, ...
    'LineWidth',2.3,'CapSize',12,'HandleVisibility','off');
set(ax,'YScale','log','XLim',[0.55 numel(x)+0.45],'XTick',x, ...
    'XTickLabel',arrayfun(@num2str,fieldSizes,'UniformOutput',false), ...
    'FontName','Arial','FontSize',26,'LineWidth',1.8,'TickDir','out','Box','off');
xlabel(ax,'Square field width (HCs)','FontSize',29);
ylabel(ax,'50-iteration end-to-end time (s)','FontSize',29);
lastRows = T.FieldHC==fieldSizes(end);
title(ax,sprintf('40 x 40 HC speedup: %.1fx', ...
    median(T.CGSeconds(lastRows)./T.DNNSeconds(lastRows))), ...
    'FontSize',28,'FontWeight','normal');
grid(ax,'on'); ax.GridAlpha=0.14;
end

function [center,error] = mean_sem(values)
values = values(:);
center = mean(values);
error = std(values,0)/sqrt(numel(values));
end
