function revise_Figure_5A_layout
% Edit the saved graphics only; do not recompute the fixed points.
outputDir = repro_paths('project/Spectral_Analysis_Report/Figures/Draft');
fig = openfig(fullfile(outputDir,'Figure_5A_layout_source.fig'),'invisible');
fig.Units = 'inches';
fig.Position = [0.3 0.3 26 26*12/22];
delete(findall(fig,'Type','legend'));
axs = findall(fig,'Type','axes');
for k=1:numel(axs)
    axs(k).Units = 'normalized';
end
positions = vertcat(axs.Position);
[~,leftIndex] = min(positions(:,1));
for k=1:numel(axs)
    if k~=leftIndex
        xlabel(axs(k),'');
        ylabel(axs(k),'');
        positions(k,1) = positions(k,1)-1.2/26;
        if any(contains(string(axs(k).Title.String),'HCnorm'))
            axs(k).Title.String = 'Condition 2';
        end
    end
end
objects = findall(fig,'-property','FontSize');
for k=1:numel(objects)
    if isprop(objects(k),'FontUnits')
        objects(k).FontUnits = 'points';
    end
    objects(k).FontSize = 45;
end
for k=1:numel(axs)
    axs(k).XLabel.FontSize = 45;
    axs(k).YLabel.FontSize = 45;
    axs(k).Title.FontSize = 45;
end
axs(leftIndex).XLabel.Position = [0.5 -0.10 0];
axs(leftIndex).YLabel.Position = [-0.13 0.5 0];
bars = findall(fig,'Type','colorbar');
for k=1:numel(bars)
    bars(k).Title.FontSize = 45;
    bars(k).Label.FontSize = 45;
    bars(k).Units = 'normalized';
    if strcmp(bars(k).Label.String,'E firing rate')
        bars(k).Position(1) = bars(k).Position(1)-1.2/26;
    end
end
drawnow;
for k=1:numel(axs)
    axs(k).Position = positions(k,:);
end
drawnow;
for k=1:numel(bars)
    if strcmp(bars(k).Label.String,'HC norm')
        bars(k).Label.Units = 'points';
        bars(k).Label.Position(1) = bars(k).Label.Position(1)-24;
    end
end
savefig(fig,fullfile(outputDir,'Figure 5A.fig'));
exportgraphics(fig,fullfile(outputDir,'Figure 5A.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding',16);
close(fig);
end
