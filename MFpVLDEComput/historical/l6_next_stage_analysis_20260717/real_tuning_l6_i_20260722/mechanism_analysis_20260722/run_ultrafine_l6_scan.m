function run_ultrafine_l6_scan()
% Resolve the narrow L6 instability around beta 0.170.

paths = mechanism_initialize();
loadedSetup = load(paths.SetupFile,'setup');
context = loadedSetup.setup.Context;
loaded = load(fullfile(paths.OutputRoot,'fine_low_branch_stability.mat'), ...
    'l6States');
fineGrid = 0.16:0.001:0.18;
state = loaded.l6States(:,find(abs(fineGrid-0.169)<1e-12,1));
betaGrid = 0.169:0.0002:0.171;
rows = cell(numel(betaGrid),1);
states = nan(numel(state),numel(betaGrid));

for index = 1:numel(betaGrid)
    beta = betaGrid(index);
    phi = @(x)mechanism_phi_variant(x,context,1+beta,1,'extended');
    fixedResult = real_tuning_fixed_point(phi,state,context.RelaxationP);
    state = fixedResult.State(:);
    states(:,index) = state;
    J = real_tuning_true_jacobian(state,context,1+beta,1);
    maxReal = mechanism_max_real(J);
    rows{index} = table(beta,maxReal,fixedResult.Iterations, ...
        fixedResult.Residual,max(state), ...
        'VariableNames',{'beta','maxRealLambda','iterations','residual', ...
        'maximumRateHz'});
    fprintf('L6 ultrafine beta=%.4f: maxRe %.9f, iter %d\n', ...
        beta,maxReal,fixedResult.Iterations);
end

ultrafineTable = vertcat(rows{:});
writetable(ultrafineTable,fullfile(paths.OutputRoot, ...
    'ultrafine_l6_stability.tsv'),'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'ultrafine_l6_stability.mat'), ...
    'ultrafineTable','states','-v7.3');
end
