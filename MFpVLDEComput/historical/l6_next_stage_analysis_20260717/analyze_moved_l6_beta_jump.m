function analyze_moved_l6_beta_jump()
% Diagnose the moved-fixed-point L6 beta jump using existing branch states.

root = fileparts(mfilename('fullpath'));
mechanismRoot = fullfile(root,'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722');
addpath(mechanismRoot,'-begin');
paths = mechanism_initialize();
loadedSetup = load(paths.SetupFile,'setup');
context = loadedSetup.setup.Context;
baseline = context.FixedPoint(:);
resultsRoot = fullfile(mechanismRoot,'results','moved_l6_beta_jump');
if ~exist(resultsRoot,'dir'); mkdir(resultsRoot); end

direct = readtable(fullfile(paths.ResultsRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
direct = sortrows(direct(direct.pathway=="L6",:),'beta');
continuation = readtable(fullfile(paths.OutputRoot,'branch_continuation.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
continuation = continuation(continuation.pathway=="L6",:);
loadedBranches = load(fullfile(paths.OutputRoot,'branch_continuation.mat'), ...
    'l6ForwardStates','l6ReverseStates');

betas = [0.275;0.300];
rows = cell(numel(betas),1);
for index = 1:numel(betas)
    beta = betas(index);
    low = local_branch_state(loadedBranches.l6ForwardStates,0:0.01:0.35, ...
        beta,context,'extended');
    high = local_branch_state(loadedBranches.l6ReverseStates,0.4:-0.01:0, ...
        beta,context,'extended');
    [lowAlpha,lowResidual] = local_metrics(low,beta,context,'extended');
    [highAlpha,highResidual] = local_metrics(high,beta,context,'extended');

    phiRaw = @(x)mechanism_phi_variant(x,context,1+beta,1,'raw');
    rawAtExtendedHighAlpha = mechanism_max_real( ...
        mechanism_jacobian_components(high,context,1+beta,1,'raw'));
    rawAtExtendedHighResidual = norm(phiRaw(high)-high)/max(norm(high),eps);
    rawLowResult = real_tuning_fixed_point(phiRaw,baseline,context.RelaxationP);
    rawHighResult = real_tuning_fixed_point(phiRaw,high,context.RelaxationP);
    rawLow = rawLowResult.State(:);
    rawHigh = rawHighResult.State(:);
    rawLowAlpha = mechanism_max_real(mechanism_jacobian_components( ...
        rawLow,context,1+beta,1,'raw'));
    rawHighAlpha = mechanism_max_real(mechanism_jacobian_components( ...
        rawHigh,context,1+beta,1,'raw'));

    rows{index} = table(beta,max(low),lowAlpha,lowResidual,max(high), ...
        highAlpha,highResidual,rawAtExtendedHighAlpha, ...
        rawAtExtendedHighResidual,max(rawLow),rawLowAlpha, ...
        rawLowResult.Residual,max(rawHigh),rawHighAlpha, ...
        rawHighResult.Residual,norm(high-low)/norm(low), ...
        norm(rawHigh-rawLow)/norm(rawLow), ...
        'VariableNames',{'beta','extendedLowMaxRate','extendedLowMaxReal', ...
        'extendedLowResidual','extendedHighMaxRate','extendedHighMaxReal', ...
        'extendedHighResidual','rawAtExtendedHighMaxReal', ...
        'rawAtExtendedHighResidual','rawLowMaxRate','rawLowMaxReal', ...
        'rawLowResidual','rawHighMaxRate','rawHighMaxReal', ...
        'rawHighResidual','extendedBranchDistance','rawBranchDistance'});
end
audit = vertcat(rows{:});
writetable(audit,fullfile(resultsRoot,'moved_l6_beta_jump_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(resultsRoot,'moved_l6_beta_jump_audit.mat'), ...
    'audit','direct','continuation','-v7.3');

local_plot(direct,continuation,audit,resultsRoot);
fprintf('%s\n',evalc('disp(audit)'));
end

function state = local_branch_state(states,grid,beta,context,responseMode)
[~,nearest] = min(abs(grid-beta));
initial = states(:,nearest);
phi = @(x)mechanism_phi_variant(x,context,1+beta,1,responseMode);
result = real_tuning_fixed_point(phi,initial,context.RelaxationP);
if ~result.Converged
    error('MovedL6:BranchRefinement','Branch refinement failed at beta %.6f.',beta);
end
state = result.State(:);
end

function [alpha,residual] = local_metrics(state,beta,context,responseMode)
phi = @(x)mechanism_phi_variant(x,context,1+beta,1,responseMode);
residual = norm(phi(state)-state)/max(norm(state),eps);
alpha = mechanism_max_real(real_tuning_true_jacobian(state,context,1+beta,1));
end

function local_plot(direct,continuation,audit,resultsRoot)
fig = figure('Color','w','Position',[100 100 1320 850],'ToolBar','none');
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on
plot(direct.beta,direct.maximumRateHz,'ko-','LineWidth',1.3, ...
    'MarkerSize',4,'DisplayName','independent runs from baseline state');
branches = unique(continuation.branch,'stable');
colors = lines(numel(branches));
for index = 1:numel(branches)
    rows = sortrows(continuation(continuation.branch==branches(index),:),'beta');
    plot(rows.beta,rows.maximumRateHz,'-','LineWidth',2, ...
        'Color',colors(index,:),'DisplayName',branches(index));
end
xline(0.275,':','HandleVisibility','off');
xline(0.300,':','HandleVisibility','off');
xlabel('\beta'); ylabel('maximum firing rate (Hz)');
title('Independent sweep versus continued fixed-point branches');
grid on; legend('Location','northwest');

nexttile; hold on
plot(direct.beta,direct.trueMaxRealLambda,'ko-','LineWidth',1.3, ...
    'MarkerSize',4,'DisplayName','independent sweep');
yline(1,'r--','DisplayName','stability boundary');
xline(0.275,':','HandleVisibility','off');
xline(0.300,':','HandleVisibility','off');
xlabel('\beta'); ylabel('max Re \lambda(D\Phi)');
title('Apparent eigenvalue drop follows branch selection');
grid on; legend('Location','best');

nexttile; hold on
plot(audit.beta,audit.extendedLowMaxRate,'o-','LineWidth',2, ...
    'DisplayName','extended moderate branch');
plot(audit.beta,audit.extendedHighMaxRate,'s-','LineWidth',2, ...
    'DisplayName','extended high branch');
plot(audit.beta,audit.rawLowMaxRate,'o--','LineWidth',1.7, ...
    'DisplayName','raw moderate branch');
plot(audit.beta,audit.rawHighMaxRate,'s--','LineWidth',1.7, ...
    'DisplayName','raw secondary branch');
xlabel('\beta'); ylabel('maximum firing rate (Hz)');
title('Coexisting roots at the jump'); grid on; legend('Location','best');

nexttile; hold on
plot(audit.beta,audit.extendedLowMaxReal,'o-','LineWidth',2, ...
    'DisplayName','extended moderate branch');
plot(audit.beta,audit.extendedHighMaxReal,'s-','LineWidth',2, ...
    'DisplayName','extended high branch');
plot(audit.beta,audit.rawLowMaxReal,'o--','LineWidth',1.7, ...
    'DisplayName','raw moderate branch');
plot(audit.beta,audit.rawHighMaxReal,'s--','LineWidth',1.7, ...
    'DisplayName','raw secondary branch');
yline(1,'r--','DisplayName','stability boundary');
xlabel('\beta'); ylabel('max Re \lambda(D\Phi)');
title('Stability of each coexisting root'); grid on; legend('Location','best');

title(layout,['Moved-fixed-point L6 beta jump: basin switch, not local ' ...
    'branch creation at \beta=0.275--0.300']);
allAxes = findall(fig,'Type','axes');
for index = 1:numel(allAxes)
    axtoolbar(allAxes(index),{});
end
exportgraphics(fig,fullfile(resultsRoot,'moved_l6_beta_jump_diagnosis.pdf'), ...
    'ContentType','vector');
savefig(fig,fullfile(resultsRoot,'moved_l6_beta_jump_diagnosis.fig'));
end
