function state = followup_branch_state(paths,pathway,beta,branch,context)
% Load the nearest saved branch state and refine it at the requested beta.

loaded = load(paths.ContinuationFile,'l6ForwardStates','l6ReverseStates', ...
    'iForwardStates','iReverseStates');
if pathway=="L6"
    if branch=="low"
        grid = 0:0.01:0.35;
        states = loaded.l6ForwardStates;
    else
        grid = 0.4:-0.01:0;
        states = loaded.l6ReverseStates;
    end
    gain6 = 1+beta;
    gainI = 1;
else
    if branch=="low"
        grid = 0:-0.005:-0.15;
        states = loaded.iForwardStates;
    else
        grid = -0.2:0.005:0;
        states = loaded.iReverseStates;
    end
    gain6 = 1;
    gainI = 1+beta;
end
[~,index] = min(abs(grid-beta));
initial = states(:,index);
phi = @(x)mechanism_phi_variant(x,context,gain6,gainI,'extended');
result = real_tuning_fixed_point(phi,initial,context.RelaxationP);
if ~result.Converged
    error('Followup:BranchState', ...
        '%s %s branch beta %+.6f failed: %s, residual %.3e.', ...
        pathway,branch,beta,result.Termination,result.Residual);
end
state = result.State(:);
end
