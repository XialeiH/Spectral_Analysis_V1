function run_continuation_stability_audit()
% Evaluate selected low/high continuation states with the exact Jacobian.

paths = mechanism_initialize();
loadedSetup = load(paths.SetupFile,'setup');
context = loadedSetup.setup.Context;
loaded = load(fullfile(paths.OutputRoot,'branch_continuation.mat'));

l6ForwardGrid = 0:0.01:0.35;
l6ReverseGrid = 0.4:-0.01:0;
iForwardGrid = 0:-0.005:-0.15;
iReverseGrid = -0.2:0.005:0;
specs = { ...
    "L6","low continuation",[0 .15 .17 .18 .20 .30 .35],l6ForwardGrid,loaded.l6ForwardStates; ...
    "L6","high continuation",[0 .15 .30 .35],l6ReverseGrid,loaded.l6ReverseStates; ...
    "Inhibition","low continuation",[0 -.05 -.075 -.08 -.10 -.125 -.145],iForwardGrid,loaded.iForwardStates; ...
    "Inhibition","high continuation",[0 -.05 -.075 -.10 -.125 -.145],iReverseGrid,loaded.iReverseStates};
rows = {};

for group = 1:size(specs,1)
    pathway = specs{group,1};
    branch = specs{group,2};
    betaList = specs{group,3};
    betaGrid = specs{group,4};
    states = specs{group,5};
    for beta = betaList
        stateIndex = find(abs(betaGrid-beta)<1e-10,1);
        state = states(:,stateIndex);
        gain6 = 1;
        gainI = 1;
        if pathway=="L6"; gain6 = 1+beta; else; gainI = 1+beta; end
        J = real_tuning_true_jacobian(state,context,gain6,gainI);
        maxReal = mechanism_max_real(J);
        stats = mechanism_domain_stats(state,context,gain6,gainI, ...
            pathway,beta,branch);
        rows{end+1,1} = table(pathway,string(branch),beta,maxReal,min(state), ...
            max(state),stats.fractionOutsideTrainingBox, ...
            stats.fractionFullExtension, ...
            'VariableNames',{'pathway','branch','beta','maxRealLambda', ...
            'minimumRateHz','maximumRateHz','fractionOutsideTrainingBox', ...
            'fractionFullExtension'}); %#ok<AGROW>
        fprintf('%s %s beta=%+.3f: maxRe=%.6f, max rate=%.3f\n', ...
            pathway,branch,beta,maxReal,max(state));
    end
end

stabilityTable = vertcat(rows{:});
writetable(stabilityTable,fullfile(paths.OutputRoot, ...
    'continuation_stability.tsv'),'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'continuation_stability.mat'),'stabilityTable','-v7.3');
end
