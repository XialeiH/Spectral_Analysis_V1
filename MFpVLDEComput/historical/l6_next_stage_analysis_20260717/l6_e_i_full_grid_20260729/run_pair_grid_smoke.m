function run_pair_grid_smoke()
% Verify baseline parity and finite moved fixed points before full submission.

setupFile=getenv('PAIR_GRID_SETUP_FILE');
outputFile=getenv('PAIR_GRID_SMOKE_OUTPUT');
if ~isfile(setupFile) || isempty(outputFile)
    error('PairGrid:SmokeEnvironment','Smoke setup and output are required.');
end
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
n=numel(baseline)/3;

oldBaseline=l6ns_phi(baseline,0,context,[1 1],1,'true');
newBaseline=l6ns_phi_three_pathway(baseline,1,1,1,context);
oldL6=l6ns_phi(baseline,1-1.1,context,[1 1],1,'true');
newL6=l6ns_phi_three_pathway(baseline,1.1,1,1,context);
baselineMapRelativeDifference=norm(newBaseline-oldBaseline)/max(norm(oldBaseline),eps);
l6MapRelativeDifference=norm(newL6-oldL6)/max(norm(oldL6),eps);

specification=[ ...
    1.0 1.0 1.0; ...
    1.1 1.0 1.0; ...
    1.0 0.9 1.0; ...
    1.0 1.1 1.0; ...
    1.0 1.0 1.1; ...
    1.1 1.1 1.0; ...
    1.1 0.9 1.0; ...
    1.0 1.1 1.1];
rowCount=size(specification,1);
converged=false(rowCount,1);
residual=nan(rowCount,1);
iterations=nan(rowCount,1);
hcNorm=nan(rowCount,1);
minimumRateHz=nan(rowCount,1);
maximumRateHz=nan(rowCount,1);
for row=1:rowCount
    gain6=specification(row,1);
    gainE=specification(row,2);
    gainI=specification(row,3);
    phi=@(state)l6ns_phi_three_pathway(state,gain6,gainE,gainI,context);
    fixed=real_tuning_fixed_point(phi,baseline,context.RelaxationP);
    converged(row)=fixed.Converged;
    residual(row)=fixed.Residual;
    iterations(row)=fixed.Iterations;
    hcNorm(row)=HC_norm_diff( ...
        fixed.State(1:n),fixed.State(n+(1:n)),fixed.State(2*n+(1:n)), ...
        baseline(1:n),baseline(n+(1:n)),baseline(2*n+(1:n)), ...
        0.3077,0.8,0.2);
    minimumRateHz(row)=min(fixed.State);
    maximumRateHz(row)=max(fixed.State);
end
smoke=table(specification(:,1),specification(:,2),specification(:,3), ...
    converged,residual,iterations,hcNorm,minimumRateHz,maximumRateHz, ...
    repmat(baselineMapRelativeDifference,rowCount,1), ...
    repmat(l6MapRelativeDifference,rowCount,1), ...
    'VariableNames',{'gain6','gainE','gainI','converged','residual', ...
    'iterations','HCnorm','minimumRateHz','maximumRateHz', ...
    'baselineMapRelativeDifference','l6MapRelativeDifference'});
writetable(smoke,outputFile,'FileType','text','Delimiter','\t');
if baselineMapRelativeDifference>1e-12 || l6MapRelativeDifference>1e-12
    error('PairGrid:Parity', ...
        'New pathway map does not reproduce the established baseline/L6 map.');
end
if ~all(converged([1:5 7]))
    error('PairGrid:SmokeConvergence', ...
        'At least one baseline or single-pathway smoke condition did not converge.');
end
disp(smoke)
end
