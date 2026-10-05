function [selectedSlope,audit] = l6ei_hc_optimal_slope(setupFile,outputRoot)
% Compute the moved-fixed-point HC-optimal betaEI-per-beta6 slope.

% For x = Phi(x,beta), u_p = (I-D_x Phi)^(-1)D_beta_p Phi.  The
% coefficient m minimizing ||u_6 + m*u_EI||_HC is used as the slope.

codeRoot = fileparts(mfilename('fullpath'));
analysisRoot = fileparts(codeRoot);
modelRoot = fullfile( ...
    '/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis', ...
    'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
runtimeRoot = fullfile(analysisRoot,'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722','runtime_h96');
addpath(modelRoot);
addpath(fullfile(modelRoot,'Utils'));
addpath(runtimeRoot,'-begin');
addpath(codeRoot,'-begin');

loaded = load(setupFile,'setup');
context = loaded.setup.Context;
x0 = context.FixedPoint(:);
j0 = sparse(loaded.setup.Pathway.JBaseline);
stateCount = numel(x0);
systemMatrix = speye(stateCount)-j0;
phi = @(state,g6,gEI)l6ns_phi_l6_ei(state,g6,gEI,context);
baselineResidual = norm(phi(x0,1,1)-x0)/max(norm(x0),eps);

steps = [1e-4;1e-5;1e-6];
stepCount = numel(steps);
slope = nan(stepCount,1);
responseCosine = nan(stepCount,1);
linearSolveRelativeResidual = nan(stepCount,1);
responses = cell(stepCount,1);
for stepIndex = 1:stepCount
    h = steps(stepIndex);
    b6 = (phi(x0,1+h,1)-phi(x0,1-h,1))/(2*h);
    bEI = (phi(x0,1,1+h)-phi(x0,1,1-h))/(2*h);
    forcing = [b6 bEI];
    response = systemMatrix\forcing;
    responses{stepIndex} = response;
    linearSolveRelativeResidual(stepIndex) = ...
        norm(systemMatrix*response-forcing,'fro')/max(norm(forcing,'fro'),eps);
    u6 = response(:,1);
    uEI = response(:,2);
    slope(stepIndex) = -local_hc_inner(uEI,u6,context)/ ...
        local_hc_inner(uEI,uEI,context);
    responseCosine(stepIndex) = local_hc_inner(u6,uEI,context)/sqrt( ...
        local_hc_inner(u6,u6,context)*local_hc_inner(uEI,uEI,context));
end

selectedIndex = 2;
selectedSlope = slope(selectedIndex);
audit = table(steps,slope,responseCosine,linearSolveRelativeResidual, ...
    repmat(baselineResidual,stepCount,1), ...
    'VariableNames',{'finiteDifferenceStep','slopeL6ItoE', ...
    'responseCosineL6ItoE','linearSolveRelativeResidual', ...
    'baselineFixedPointRelativeResidual'});
writetable(audit,fullfile(outputRoot, ...
    'theoretical_hc_optimal_fixedpoint_slope_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
selectedResponse = responses{selectedIndex};
save(fullfile(outputRoot, ...
    'theoretical_hc_optimal_fixedpoint_slope_audit.mat'), ...
    'audit','selectedSlope','selectedResponse','-v7.3');
fprintf(['L6/I-to-E moved-fixed-point HC slope %.12g; cosine %.6f; ' ...
    'baseline residual %.3e.\n'],selectedSlope, ...
    responseCosine(selectedIndex),baselineResidual);
end

function value = local_hc_inner(u,v,context)
n = numel(u)/3;
uE = (1-context.CWeight)*u(1:n)+context.CWeight*u(n+(1:n));
vE = (1-context.CWeight)*v(1:n)+context.CWeight*v(n+(1:n));
uI = u(2*n+(1:n));
vI = v(2*n+(1:n));
value = real(mean(0.8*uE.*vE+0.2*uI.*vI));
end
