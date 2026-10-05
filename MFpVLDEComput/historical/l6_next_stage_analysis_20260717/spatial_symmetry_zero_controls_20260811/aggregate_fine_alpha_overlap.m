function summary = aggregate_fine_alpha_overlap(outputRoot)
% Aggregate the fine-alpha overlap diagnostic and plot threshold counts.

if nargin < 1 || isempty(outputRoot)
    outputRoot = getenv('SPATIAL_CONTROL_OUTPUT');
end
files = dir(fullfile(outputRoot,'fine_alpha_*.tsv'));
files = files(~strcmp({files.name},'fine_alpha_overlap_summary.tsv'));
if numel(files) ~= 17
    error('FineAlpha:Missing','Expected 17 task summaries, found %d.',numel(files));
end
parts = cell(numel(files),1);
for index = 1:numel(files)
    parts{index} = readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t');
end
summary = sortrows(vertcat(parts{:}),'alpha');
summary.schurMedianNearestError = NaN(height(summary),1);
summary.schurP95NearestError = NaN(height(summary),1);
summary.generalizedMedianNearestError = NaN(height(summary),1);
summary.generalizedP95NearestError = NaN(height(summary),1);
summary.generalizedMedianNearestRelativeError = NaN(height(summary),1);
summary.generalizedP95NearestRelativeError = NaN(height(summary),1);
for row = 1:height(summary)
    tag = strrep(sprintf('%.4f',summary.alpha(row)),'.','p');
    loaded = load(fullfile(outputRoot,['fine_alpha_' tag '.mat']), ...
        'fullEigenvalues','schurEigenvalues','generalizedSchurEigenvalues');
    [schurError,~] = local_nearest_error( ...
        loaded.schurEigenvalues,loaded.fullEigenvalues);
    [generalizedError,generalizedRelativeError] = local_nearest_error( ...
        loaded.generalizedSchurEigenvalues,loaded.fullEigenvalues);
    summary.schurMedianNearestError(row) = median(schurError);
    summary.schurP95NearestError(row) = prctile(schurError,95);
    summary.generalizedMedianNearestError(row) = median(generalizedError);
    summary.generalizedP95NearestError(row) = prctile(generalizedError,95);
    summary.generalizedMedianNearestRelativeError(row) = ...
        median(generalizedRelativeError);
    summary.generalizedP95NearestRelativeError(row) = ...
        prctile(generalizedRelativeError,95);
end
writetable(summary,fullfile(outputRoot,'fine_alpha_overlap_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Color','w','Position',[100 100 1380 980]);
layout = tiledlayout(fig,2,2,'TileSpacing','loose','Padding','compact');
local_panel(nexttile(layout),summary,'fullSigmaLE0p01','fullSigmaLE0p025', ...
    'fullSigmaLE0p05','Full Jacobian: small singular values');
local_panel(nexttile(layout),summary,'fullEigLE0p01','fullEigLE0p025', ...
    'fullEigLE0p05','Full Jacobian: near-zero eigenvalues');
local_panel(nexttile(layout),summary,'decoupledEigLE0p01', ...
    'decoupledEigLE0p025','decoupledEigLE0p05', ...
    'E-I feedback removed: near-zero eigenvalues');
local_panel(nexttile(layout),summary,'schurEigLE0p01', ...
    'schurEigLE0p025','schurEigLE0p05', ...
    'E-I Schur cancellation matrix H');
title(layout,'Spatial smoothing threshold crossings and E-I overlap');
exportgraphics(fig,fullfile(outputRoot,'06_fine_alpha_gaussian_ei_overlap.pdf'), ...
    'ContentType','vector');
savefig(fig,fullfile(outputRoot,'06_fine_alpha_gaussian_ei_overlap.fig'));
end

function [nearestError,nearestRelativeError] = local_nearest_error(approximate,direct)
nearestError = zeros(numel(approximate),1);
nearestRelativeError = zeros(numel(approximate),1);
for index = 1:numel(approximate)
    [nearestError(index),nearestIndex] = min(abs(direct-approximate(index)));
    nearestRelativeError(index) = nearestError(index)/ ...
        max(abs(direct(nearestIndex)),1e-3);
end
end

function local_panel(ax,summary,field1,field2,field3,panelTitle)
plot(ax,summary.alpha,summary.(field1),'-o','LineWidth',1.5,'MarkerSize',4); hold(ax,'on');
plot(ax,summary.alpha,summary.(field2),'-s','LineWidth',1.5,'MarkerSize',4);
plot(ax,summary.alpha,summary.(field3),'-^','LineWidth',1.5,'MarkerSize',4);
xlabel(ax,'joint spatial-flattening strength \alpha');
ylabel(ax,'count'); title(ax,panelTitle); grid(ax,'on'); box(ax,'on');
ax.FontSize = 10;
legend(ax,{'threshold 0.01','threshold 0.025','threshold 0.05'}, ...
    'Location','best');
end
