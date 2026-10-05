function summary = aggregate_spatial_smoothing_rank_controls(outputRoot)
files = dir(fullfile(outputRoot,'summary_*.tsv'));
if numel(files)~=7
    error('SpatialControl:Missing','Expected 7 summaries, found %d.',numel(files));
end
parts = cell(numel(files),1);
for index=1:numel(files)
    parts{index}=readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t');
end
summary=vertcat(parts{:});
orderNames=["baseline";"joint_alpha_0p25";"joint_alpha_0p50"; ...
    "joint_alpha_0p75";"L6_delta";"L4_delta";"L4_L6_delta"];
[~,order]=ismember(string(summary.control),orderNames);
[~,sortIndex]=sort(order);
summary=summary(sortIndex,:);
writetable(summary,fullfile(outputRoot,'spatial_smoothing_rank_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig=figure('Visible','off','Color','w','Position',[50 50 1500 760]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);
bar(ax,[summary.countAbsEigLE0p05 summary.countSingularLE0p05]);
grid(ax,'on'); box(ax,'on'); ylabel(ax,'Count at threshold 0.05');
xticks(ax,1:height(summary)); xticklabels(ax,replace(string(summary.control),'_',' '));
xtickangle(ax,25); legend(ax,{'|lambda| <= 0.05','sigma <= 0.05'}, ...
    'Location','northoutside','Orientation','horizontal');
title(ax,'Near-zero eigenvalues and small singular values');
ax=nexttile(layout);
yyaxis(ax,'left'); plot(ax,1:height(summary),summary.stableRank,'o-','LineWidth',1.8);
ylabel(ax,'Stable rank');
yyaxis(ax,'right'); plot(ax,1:height(summary),summary.entropyEffectiveRank, ...
    's-','LineWidth',1.8); ylabel(ax,'Entropy effective rank');
xticks(ax,1:height(summary)); xticklabels(ax,replace(string(summary.control),'_',' '));
xtickangle(ax,25); grid(ax,'on'); box(ax,'on');
title(ax,'Effective rank after flattening spatial smoothing');
title(layout,{'Spatial-smoothing control at fixed point and fixed local gains', ...
    'All controls preserve block row sums and spatial symmetries'});
exportgraphics(fig,fullfile(outputRoot,'05_spatial_smoothing_effective_rank_control.pdf'), ...
    'ContentType','vector');
close(fig);
end
