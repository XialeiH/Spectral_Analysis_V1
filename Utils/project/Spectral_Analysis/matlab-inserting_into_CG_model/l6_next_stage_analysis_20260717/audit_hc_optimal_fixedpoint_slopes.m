function slopes = audit_hc_optimal_fixedpoint_slopes(setupFile, outputRoot)
% Compute first-order pathway compensation slopes in the canonical HC norm.

% For x = Phi(x,beta), the fixed-point sensitivity is
% u_p = (I-D_x Phi)^(-1) D_{beta_p} Phi.  For beta_q=m*beta_p,
% m minimizes ||u_p+m*u_q||_HC.

if nargin < 1 || isempty(setupFile)
    codeRoot = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');
    setupFile = fullfile(codeRoot, 'results_global_bifurcation_20260722', ...
        'global_bifurcation_setup.mat');
else
    codeRoot = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile( ...
        [repro_paths('project') '/Spectral_Analysis'], ...
        'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
        'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors', ...
        'L6 and Inhibition');
end

modelRoot = fullfile( ...
    [repro_paths('project') '/Spectral_Analysis'], ...
    'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
runtimeRoot = fullfile(codeRoot,'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722','runtime_h96');
repro_addpath(modelRoot);
repro_addpath(fullfile(modelRoot,'Utils'));
repro_addpath(runtimeRoot,'-begin');
repro_addpath(fullfile(codeRoot,'l6_e_i_full_grid_20260729'),'-begin');
loaded = load(setupFile, 'setup');
setup = loaded.setup;
context = setup.Context;
x0 = context.FixedPoint(:);
j0 = sparse(setup.Pathway.JBaseline);
stateCount = numel(x0);

phi = @(state,g6,gE,gI)l6ns_phi_three_pathway( ...
    state,g6,gE,gI,context);
baselineResidual = norm(phi(x0,1,1,1)-x0) / max(norm(x0),eps);

steps = [1e-4; 1e-5; 1e-6];
stepCount = numel(steps);
slope6I = zeros(stepCount,1);
slope6E = zeros(stepCount,1);
slopeEI = zeros(stepCount,1);
cosine6I = zeros(stepCount,1);
cosine6E = zeros(stepCount,1);
cosineEI = zeros(stepCount,1);
solveResidual = zeros(stepCount,1);
responses = cell(stepCount,1);

systemMatrix = speye(stateCount)-j0;
for stepIndex = 1:stepCount
    h = steps(stepIndex);
    b6 = (phi(x0,1+h,1,1)-phi(x0,1-h,1,1))/(2*h);
    bE = (phi(x0,1,1+h,1)-phi(x0,1,1-h,1))/(2*h);
    bI = (phi(x0,1,1,1+h)-phi(x0,1,1,1-h))/(2*h);
    forcing = [b6 bE bI];
    response = systemMatrix\forcing;
    responses{stepIndex} = response;
    solveResidual(stepIndex) = norm(systemMatrix*response-forcing,'fro') / ...
        max(norm(forcing,'fro'),eps);

    u6 = response(:,1);
    uE = response(:,2);
    uI = response(:,3);
    slope6I(stepIndex) = -local_hc_inner(uI,u6,context) / ...
        local_hc_inner(uI,uI,context);
    slope6E(stepIndex) = -local_hc_inner(uE,u6,context) / ...
        local_hc_inner(uE,uE,context);
    slopeEI(stepIndex) = -local_hc_inner(uI,uE,context) / ...
        local_hc_inner(uI,uI,context);
    cosine6I(stepIndex) = local_hc_cosine(u6,uI,context);
    cosine6E(stepIndex) = local_hc_cosine(u6,uE,context);
    cosineEI(stepIndex) = local_hc_cosine(uE,uI,context);
end

selectedIndex = 2;
selectedResponse = responses{selectedIndex};
slopes = struct( ...
    'L6Inhibition',slope6I(selectedIndex), ...
    'L6Excitation',slope6E(selectedIndex), ...
    'ExcitationInhibition',slopeEI(selectedIndex), ...
    'CosineL6Inhibition',cosine6I(selectedIndex), ...
    'CosineL6Excitation',cosine6E(selectedIndex), ...
    'CosineExcitationInhibition',cosineEI(selectedIndex), ...
    'FiniteDifferenceStep',steps(selectedIndex));

audit = table(steps,slope6I,slope6E,slopeEI,cosine6I,cosine6E, ...
    cosineEI,solveResidual,repmat(baselineResidual,stepCount,1), ...
    'VariableNames',{'finiteDifferenceStep','slopeL6Inhibition', ...
    'slopeL6Excitation','slopeExcitationInhibition', ...
    'cosineL6Inhibition','cosineL6Excitation', ...
    'cosineExcitationInhibition','linearSolveRelativeResidual', ...
    'baselineFixedPointRelativeResidual'});
writetable(audit, fullfile(outputRoot, ...
    'theoretical_hc_optimal_fixedpoint_slope_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot, ...
    'theoretical_hc_optimal_fixedpoint_slope_audit.mat'), ...
    'audit','slopes','selectedResponse','-v7.3');

fprintf(['HC-optimal fixed-point slopes: L6-I %.12g, L6-E %.12g, ' ...
    'E-I %.12g; response cosines %.6f, %.6f, %.6f.\n'], ...
    slopes.L6Inhibition,slopes.L6Excitation, ...
    slopes.ExcitationInhibition,slopes.CosineL6Inhibition, ...
    slopes.CosineL6Excitation,slopes.CosineExcitationInhibition);
end

function value = local_hc_inner(u,v,context)
n = numel(u)/3;
uE = (1-context.CWeight)*u(1:n)+context.CWeight*u(n+(1:n));
vE = (1-context.CWeight)*v(1:n)+context.CWeight*v(n+(1:n));
uI = u(2*n+(1:n));
vI = v(2*n+(1:n));
value = real(mean(0.8*uE.*vE+0.2*uI.*vI));
end

function value = local_hc_cosine(u,v,context)
value = local_hc_inner(u,v,context) / sqrt( ...
    local_hc_inner(u,u,context)*local_hc_inner(v,v,context));
end
