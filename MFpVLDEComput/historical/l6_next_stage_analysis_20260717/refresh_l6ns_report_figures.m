function refresh_l6ns_report_figures(resultRoot)
% Refresh report-only labels and ordering without recomputing analysis data.

if nargin < 1 || isempty(resultRoot)
    resultRoot = fullfile(fileparts(mfilename('fullpath')), ...
        'results_two_pathway_final');
end

trueRoot = fullfile(resultRoot,'true_parameter_validation');
trueStem = fullfile(trueRoot,'figure6_true_vs_fpp_validation');
fig = openfig([trueStem '.fig'],'invisible');
axesHandles = local_sorted_axes(fig);
for index = 1:numel(axesHandles)
    lines = findobj(axesHandles(index),'Type','line');
    for lineIndex = 1:numel(lines)
        x = lines(lineIndex).XData;
        y = lines(lineIndex).YData;
        [x,order] = sort(x);
        lines(lineIndex).XData = x;
        lines(lineIndex).YData = y(order);
    end
end
local_resave(fig,trueStem);

planRoot = fullfile(resultRoot,'two_pathway_plan');
fateStem = fullfile(planRoot,'figure6_inhibition_first_boundary_fate');
fig = openfig([fateStem '.fig'],'invisible');
ax = findobj(fig,'Type','axes');
lines = flipud(findobj(ax,'Type','line'));
labels = {'below boundary (-)','below boundary (+)', ...
    'above boundary (-)','above boundary (+)'};
for index = 1:min(numel(lines),numel(labels))
    lines(index).DisplayName = labels{index};
end
oldLegend = findobj(fig,'Type','legend');
delete(oldLegend);
legend(ax,lines(1:min(4,numel(lines))),labels(1:min(4,numel(lines))), ...
    'Location','best');
boundary = readtable(fullfile(planRoot,'inhibition_first_boundary.tsv'), ...
    'FileType','text','Delimiter','\t');
title(ax,sprintf('First inhibition boundary: gamma_I^c = %.6f', ...
    boundary.gammaICritical(1)));
local_resave(fig,fateStem);

local_label_map_rows(fullfile(planRoot, ...
    'figure4_l6_inhibition_left_right_maps'), ...
    {'L6 right','L6 left','I right','I left'});
local_label_map_rows(fullfile(planRoot, ...
    'figure3_matched_singular_input_output_maps'), ...
    {'L6 input','L6 output','I input','I output'});
end

function local_label_map_rows(stem,labels)
fig = openfig([stem '.fig'],'invisible');
axesHandles = local_sorted_axes(fig);
delete(findall(fig,'Tag','L6NSRowLabel'));
for row = 1:min(numel(labels),floor(numel(axesHandles)/4))
    ax = axesHandles(4*(row-1)+1);
    text(ax,-0.12,0.5,labels{row},'Units','normalized', ...
        'Rotation',90,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontWeight','bold', ...
        'Tag','L6NSRowLabel');
end
local_resave(fig,stem);
end

function axesHandles = local_sorted_axes(fig)
axesHandles = findobj(fig,'Type','axes');
positions = vertcat(axesHandles.Position);
[~,order] = sortrows(positions,[-2 1]);
axesHandles = axesHandles(order);
end

function local_resave(fig,stem)
savefig(fig,[stem '.fig']);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
exportgraphics(fig,[stem '.png'],'Resolution',220);
close(fig);
end
