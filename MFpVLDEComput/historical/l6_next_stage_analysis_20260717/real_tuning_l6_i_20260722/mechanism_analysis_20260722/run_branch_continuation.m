function run_branch_continuation()
% Follow low- and high-rate fixed-point branches in both beta directions.

paths = mechanism_initialize();
loaded = load(paths.SetupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
summary = readtable(fullfile(paths.ResultsRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');

l6High = local_saved_state(paths.ResultsRoot,summary,"L6",0.4);
iHigh = local_saved_state(paths.ResultsRoot,summary,"Inhibition",-0.2);

[l6ForwardTable,l6ForwardStates] = local_follow( ...
    "L6","low branch, increasing beta",0:0.01:0.35,baseline,context);
[l6ReverseTable,l6ReverseStates] = local_follow( ...
    "L6","high branch, decreasing beta",0.4:-0.01:0,l6High,context);
[iForwardTable,iForwardStates] = local_follow( ...
    "Inhibition","low branch, decreasing beta",0:-0.005:-0.15,baseline,context);
[iReverseTable,iReverseStates] = local_follow( ...
    "Inhibition","high branch, increasing beta",-0.2:0.005:0,iHigh,context);

continuationTable = [l6ForwardTable;l6ReverseTable;iForwardTable;iReverseTable];
writetable(continuationTable,fullfile(paths.OutputRoot,'branch_continuation.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'branch_continuation.mat'), ...
    'continuationTable','l6ForwardStates','l6ReverseStates', ...
    'iForwardStates','iReverseStates','-v7.3');
local_plot(continuationTable,paths.FigureRoot);
end

function [rows,states] = local_follow(pathway,branch,betaGrid,initialState,context)
states = nan(numel(initialState),numel(betaGrid));
rowCells = cell(numel(betaGrid),1);
state = initialState(:);
baseline = context.FixedPoint(:);
for index = 1:numel(betaGrid)
    beta = betaGrid(index);
    gain6 = 1;
    gainI = 1;
    if pathway=="L6"; gain6 = 1+beta; else; gainI = 1+beta; end
    phi = @(x)mechanism_phi_variant(x,context,gain6,gainI,'extended');
    result = real_tuning_fixed_point(phi,state,context.RelaxationP);
    state = result.State(:);
    states(:,index) = state;
    relativeShift = norm(state-baseline)/max(norm(baseline),eps);
    domain = mechanism_domain_stats(state,context,gain6,gainI, ...
        pathway,beta,branch);
    rowCells{index} = table(pathway,string(branch),beta,gain6,gainI, ...
        result.Converged,string(result.Termination),result.Iterations, ...
        result.Residual,min(state),max(state),relativeShift, ...
        domain.fractionOutsideTrainingBox,domain.fractionFullExtension, ...
        domain.rawVsExtendedRelativeResponseDifference, ...
        'VariableNames',{'pathway','branch','beta','gain6','gainI', ...
        'converged','termination','iterations','residual','minimumRateHz', ...
        'maximumRateHz','relativeBaselineShift','fractionOutsideTrainingBox', ...
        'fractionFullExtension','rawVsExtendedRelativeResponseDifference'});
    fprintf('%s | %s | beta=%+.3f: max=%.3f, outside=%.3f, iter=%d, residual=%.2e\n', ...
        pathway,branch,beta,max(state),domain.fractionOutsideTrainingBox, ...
        result.Iterations,result.Residual);
    if ~result.Converged
        warning('Mechanism:Continuation','Continuation did not converge at beta %+.3f.',beta);
    end
end
rows = vertcat(rowCells{:});
end

function state = local_saved_state(resultsRoot,summary,pathway,beta)
mask = summary.pathway==pathway & abs(summary.beta-beta)<1e-10;
if nnz(mask)~=1
    error('Mechanism:CaseLookup','Could not resolve %s beta %+.3f.',pathway,beta);
end
taskId = summary.taskId(mask);
match = dir(fullfile(resultsRoot,sprintf('%03d_*',taskId),'case_result.mat'));
loaded = load(fullfile(match.folder,match.name),'moved');
state = loaded.moved(:);
end

function local_plot(tableData,figureRoot)
pathways = ["L6","Inhibition"];
figure('Color','w','Position',[100 100 1260 500]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for index = 1:2
    nexttile; hold on
    subset = tableData(tableData.pathway==pathways(index),:);
    branches = unique(subset.branch,'stable');
    for branchIndex = 1:numel(branches)
        rows = subset(subset.branch==branches(branchIndex),:);
        [~,order] = sort(rows.beta);
        rows = rows(order,:);
        plot(rows.beta,rows.maximumRateHz,'-o','LineWidth',1.5, ...
            'DisplayName',branches(branchIndex));
    end
    xlabel('\beta'); ylabel('maximum fixed-point firing rate (Hz)');
    title(sprintf('%s forward/reverse continuation',pathways(index)));
    grid on; legend('Location','best');
end
title(layout,'Branch selection and hysteresis under true pathway tuning');
exportgraphics(gcf,fullfile(figureRoot,'03_branch_continuation_hysteresis.pdf'), ...
    'ContentType','vector');
end
