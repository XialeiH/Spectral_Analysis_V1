function smoke_all_spatial_controls_reequilibrated(setupFile)
% Validate controlled nonlinear maps and exact baseline Jacobian parity.

if nargin<1 || isempty(setupFile); setupFile=getenv('SPATIAL_CONTROL_SETUP'); end
loaded=load(setupFile,'setup');
setup=loaded.setup;
context=setup.Context;
fixedPoint=context.FixedPoint(:);
controls=["baseline";"shared_source_permutation"; ...
    "joint_alpha_0p10";"joint_alpha_0p25";"joint_alpha_0p50"; ...
    "remove_L6_smoothing";"remove_L4_smoothing"; ...
    "remove_L4_and_L6_smoothing";"L4_left_half";"L4_triangle"];

fprintf('Iteration helper: %s\n',which( ...
    'LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle'));
for index=1:numel(controls)
    [operators,metadata]=l6ns_full_spatial_control_operators(context,controls(index));
    response=l6ns_controlled_phi(fixedPoint,context,operators);
    assert(all(isfinite(response)),'%s produced nonfinite response.',controls(index));
    fprintf('%-32s initial residual %.6e, row errors %.3e %.3e\n', ...
        controls(index),norm(response-fixedPoint)/max(norm(fixedPoint),eps), ...
        metadata.MaximumRowSumRelativeError,metadata.L6RowSumRelativeError);
end

[baselineOperators,~]=l6ns_full_spatial_control_operators(context,"baseline");
J=l6ns_controlled_jacobian(setup,baselineOperators,fixedPoint);
jacobianRelativeError=norm(J-setup.Pathway.JBaseline,'fro')/ ...
    max(norm(setup.Pathway.JBaseline,'fro'),eps);
baselineResidual=norm(l6ns_controlled_phi(fixedPoint,context,baselineOperators)-fixedPoint)/ ...
    max(norm(fixedPoint),eps);
fprintf('Baseline residual %.12e\n',baselineResidual);
fprintf('Baseline Jacobian relative error %.12e\n',jacobianRelativeError);
assert(baselineResidual<1e-5, ...
    'Baseline response map does not reproduce the saved equilibrium cache.');
assert(jacobianRelativeError<1e-8,'Baseline controlled Jacobian does not reproduce JBaseline.');
fprintf('All spatial-control smoke checks passed.\n');
end
