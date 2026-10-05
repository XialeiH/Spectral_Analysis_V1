function aggregate_real_moved_iteration_results()
% Aggregate real moved-fixed-point regular-iteration tasks.

setupFile = getenv('REAL_MOVED_SETUP_FILE');
outputRoot = getenv('REAL_MOVED_OUTPUT_ROOT');
if isempty(setupFile) || ~isfile(setupFile)
    error('RealMoved:SetupFile', 'REAL_MOVED_SETUP_FILE is missing.');
end
if isempty(outputRoot) || ~isfolder(outputRoot)
    error('RealMoved:OutputRoot', 'REAL_MOVED_OUTPUT_ROOT is missing.');
end
loaded = load(setupFile, 'setup');
setup = loaded.setup;

caseFiles = dir(fullfile(outputRoot, '*', 'case_summary.tsv'));
snapshotFiles = dir(fullfile(outputRoot, '*', 'snapshot_summary.tsv'));
if numel(caseFiles) ~= 10 || numel(snapshotFiles) ~= 10
    error('RealMoved:Incomplete', ...
        'Expected 10 cases, found %d case and %d snapshot files.', ...
        numel(caseFiles), numel(snapshotFiles));
end
caseTables = cell(numel(caseFiles),1);
for index = 1:numel(caseFiles)
    caseTables{index} = readtable(fullfile(caseFiles(index).folder, ...
        caseFiles(index).name), 'FileType','text','Delimiter','\t', ...
        'TextType','string');
end
caseSummary = sortrows(vertcat(caseTables{:}),'taskId');
snapshotTables = cell(numel(snapshotFiles),1);
for index = 1:numel(snapshotFiles)
    snapshotTables{index} = readtable(fullfile(snapshotFiles(index).folder, ...
        snapshotFiles(index).name), 'FileType','text','Delimiter','\t', ...
        'TextType','string');
end
snapshotSummary = sortrows(vertcat(snapshotTables{:}), ...
    {'taskId','snapshotIndex'});
writetable(caseSummary,fullfile(outputRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(snapshotSummary,fullfile(outputRoot,'all_snapshot_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

local_plot_baseline(setup.Context.FixedPoint(:), ...
    double(setup.Context.MapSize(:).'), setup.Context.CWeight, setup, outputRoot);
local_plot_hcnorm_summary(caseSummary,snapshotSummary,setup,outputRoot);
local_write_manifest(outputRoot,caseSummary,snapshotSummary);

fprintf('Aggregated %d real moved-fixed-point iteration cases under %s.\n', ...
    height(caseSummary),outputRoot);
disp(caseSummary(:,{'pathway','beta','pathwayGain','fixedPointIterations', ...
    'fixedPointResidual','maxRealLambda','movedHCnormFromOriginalBaseline', ...
    'iterationReplayRelativeError'}));
end

function local_plot_baseline(fixed,mapSize,wC,setup,outputRoot)
maps = local_state_maps(fixed,mapSize,wC);
names = {'S','C','I','E=(1-0.3077)S+0.3077C'};
fig = figure('Visible','off','Color','w','Position',[100 100 1200 1000]);
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for index = 1:4
    ax = nexttile(layout,index);
    imagesc(ax,maps{index});
    axis(ax,'image'); axis(ax,'off');
    title(ax,names{index},'FontWeight','bold');
    colorbar(ax,'eastoutside');
end
colormap(fig,parula(256));
title(layout,sprintf(['Original baseline fixed point before real pathway tuning | ' ...
    'angle %.1f deg, contrast %g'], ...
    setup.Config.Angle,setup.Config.Contrast),'FontWeight','bold');
exportgraphics(fig,fullfile(outputRoot, ...
    '00_original_baseline_fixed_point_S_C_I_E.pdf'), ...
    'ContentType','image','Resolution',220);
close(fig);
end

function local_plot_hcnorm_summary(caseSummary,snapshotSummary,setup,outputRoot)
fig = figure('Visible','off','Color','w','Position',[100 100 1400 660]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
pathways = ["L6","Inhibition"];
colors = lines(5);
for pathwayIndex = 1:2
    ax = nexttile(layout,pathwayIndex);
    hold(ax,'on');
    cases = caseSummary(caseSummary.pathway==pathways(pathwayIndex),:);
    cases = sortrows(cases,'taskId');
    for caseIndex = 1:height(cases)
        rows = snapshotSummary(snapshotSummary.taskId==cases.taskId(caseIndex),:);
        plot(ax,rows.iteration,rows.HCnormFromOriginalBaseline,'-o', ...
            'LineWidth',1.7,'MarkerSize',4,'Color',colors(caseIndex,:), ...
            'DisplayName',sprintf(['beta=%+.3f, gain=%.3f, ' ...
            'max Re lambda=%.3f'], cases.beta(caseIndex), ...
            cases.pathwayGain(caseIndex),cases.maxRealLambda(caseIndex)));
    end
    set(ax,'YScale','log');
    grid(ax,'on'); box(ax,'on');
    xlabel(ax,'Regular iteration');
    ylabel(ax,'HCnorm from original baseline fixed point');
    title(ax,sprintf('%s real pathway iteration',pathways(pathwayIndex)));
    legend(ax,'Location','best','FontSize',8);
end
title(layout,sprintf(['Real moved-fixed-point regular iteration | ' ...
    'relaxation p=%.2f | every endpoint is the recomputed fixed point'], ...
    setup.Context.RelaxationP),'FontWeight','bold');
exportgraphics(fig,fullfile(outputRoot, ...
    '11_HCnorm_trajectories_real_L6_and_Inhibition_iteration.pdf'), ...
    'ContentType','vector');
close(fig);
end

function maps = local_state_maps(state,mapSize,wC)
n = prod(mapSize);
maps = cell(4,1);
maps{1} = reshape(state(1:n),mapSize);
maps{2} = reshape(state(n+(1:n)),mapSize);
maps{3} = reshape(state(2*n+(1:n)),mapSize);
maps{4} = (1-wC)*maps{1}+wC*maps{2};
end

function local_write_manifest(outputRoot,caseSummary,snapshotSummary)
pdfFiles = dir(fullfile(outputRoot,'**','*.pdf'));
paths = strings(numel(pdfFiles),1);
for index = 1:numel(pdfFiles)
    paths(index) = string(fullfile(pdfFiles(index).folder,pdfFiles(index).name));
end
manifest = table(paths,'VariableNames',{'pdfFile'});
writetable(manifest,fullfile(outputRoot,'figure_manifest.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'real_moved_iteration_summary.mat'), ...
    'caseSummary','snapshotSummary','manifest','-v7.3');
end
