function summary = gke_plot_static_rank_bars(dataRoot,pureJitterRoot,outputFile,summaryFile)
% Compare exact singular-spectrum ranks of equilibrium Jacobians.

catalog = gke_condition_catalog();
shortLabels = ["Baseline","E/I ranges exchanged", ...
    "Triangle","Square","Hexagon","Left half-kernel", ...
    "Joint alpha=0.10","Joint alpha=0.25","Joint alpha=0.50", ...
    "No L6 smoothing","No L4 smoothing","No L4/L6 smoothing", ...
    "Jitter 1 px","Pure jitter 2 px","Sinusoidal warp 2 px", ...
    "25% source permutation","Kernel noise 0.10","Kernel noise 0.30", ...
    "Strength noise 0.05","Strength noise 0.15"];

entropyEffectiveRank = zeros(height(catalog),1);
participationRatio = zeros(height(catalog),1);
for catalogIndex = 1:height(catalog)
    inputFile = fullfile(dataRoot,catalog.name(catalogIndex)+".mat");
    if catalog.name(catalogIndex) == "translation_jitter_2px_mix0p50"
        inputFile = fullfile(pureJitterRoot,catalog.name(catalogIndex)+".mat");
    end
    loaded = load(inputFile,'singularValues');
    singularValues = double(loaded.singularValues(:));
    singularValues = singularValues(singularValues >= 0);

    nuclearWeights = singularValues/sum(singularValues);
    entropyEffectiveRank(catalogIndex) = exp( ...
        -sum(nuclearWeights.*log(max(nuclearWeights,eps))));

    energyWeights = singularValues.^2/sum(singularValues.^2);
    participationRatio(catalogIndex) = 1/sum(energyWeights.^2);
end

summary = table(catalog.id,catalog.name,catalog.group,shortLabels(:), ...
    entropyEffectiveRank,participationRatio, ...
    'VariableNames',{'id','name','group','label', ...
    'entropyEffectiveRank','participationRatio'});
writetable(summary,summaryFile,'FileType','text','Delimiter','\t');

colors = turbo(height(catalog));
colors(1,:) = [0.10 0.10 0.10];
figureHandle = figure('Color','w','Position',[20 20 2200 1900]);
axes('Parent',figureHandle,'Position',[0.03 0.55 0.94 0.37]);
local_rank_barh(entropyEffectiveRank,colors,shortLabels, ...
    'Entropy effective rank', ...
    'Entropy rank: exp[-sum(p_i log(p_i))],  p_i = sigma_i / sum(sigma)');

axes('Parent',figureHandle,'Position',[0.03 0.07 0.94 0.37]);
local_rank_barh(participationRatio,colors,shortLabels, ...
    'Singular-spectrum participation ratio', ...
    'Participation ratio: 1 / sum(q_i^2),  q_i = sigma_i^2 / sum(sigma^2)');

annotation(figureHandle,'textbox',[0.03 0.965 0.94 0.025], ...
    'String','9.20 Static effective ranks of equilibrium Jacobians', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontWeight','bold','FontSize',19);
set(findall(figureHandle,'Type','axes'),'FontSize',11,'LineWidth',1.0);
exportgraphics(figureHandle,outputFile,'ContentType','vector');
close(figureHandle);
end

function local_rank_barh(values,colors,labels,xLabelText,titleText)
barHandle = barh(values,'FaceColor','flat','EdgeColor',[0.20 0.20 0.20], ...
    'LineWidth',0.7);
barHandle.CData = colors;
hold on
xline(values(1),'--','Color',[0.15 0.15 0.15],'LineWidth',1.5);
grid on
labelWidth = 0.18*max(values);
xlim([-labelWidth 1.17*max(values)]);
ylim([0.3 numel(values)+0.7]);
yticks(1:numel(labels));
yticklabels(repmat({''},numel(labels),1));
set(gca,'YDir','reverse');
xlabel(xLabelText,'FontSize',14);
title(titleText,'Interpreter','none','FontSize',13);
for index = 1:numel(values)
    text(-0.97*labelWidth,index,labels(index), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontSize',10,'Color',[0.10 0.10 0.10]);
    text(values(index)+0.012*max(values),index,sprintf('%.1f',values(index)), ...
        'HorizontalAlignment','left','VerticalAlignment','middle', ...
        'FontSize',10);
end
end
