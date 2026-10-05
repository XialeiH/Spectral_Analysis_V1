function result = ode_tuning_fixed_point(phi, initialState)
% Continue the model ODE until it reaches a fixed point or a hard limit.

state = initialState(:);
maximumIntrinsicTime = 1000;
chunkIntrinsicTime = 50;
tolerance = 1e-9;
rateLimitHz = 500;
intrinsicTime = 0;
acceptedSteps = 0;
termination = "maximum_time";

options = odeset('RelTol',2e-6,'AbsTol',1e-8,'MaxStep',5, ...
    'Refine',1,'Events',@local_rate_event);
residual = local_residual(state);
while residual>tolerance && intrinsicTime<maximumIntrinsicTime
    duration = min(chunkIntrinsicTime,maximumIntrinsicTime-intrinsicTime);
    solution = ode45(@(~,value)phi(value)-value,[0 duration],state,options);
    state = solution.y(:,end);
    intrinsicTime = intrinsicTime+solution.x(end);
    acceptedSteps = acceptedSteps+numel(solution.x)-1;
    if any(~isfinite(state))
        termination = "nonfinite";
        break
    end
    if max(abs(state))>=rateLimitHz
        termination = "rate_limit";
        break
    end
    residual = local_residual(state);
    if solution.x(end)<duration-1e-10*max(1,duration)
        termination = "rate_limit";
        break
    end
end
if residual<=tolerance
    termination = "converged";
end

result = struct('State',state,'Residual',residual, ...
    'IntrinsicTime',intrinsicTime,'AcceptedSteps',acceptedSteps, ...
    'Converged',residual<=tolerance,'Termination',termination, ...
    'RateLimitHz',rateLimitHz);

    function value = local_residual(candidate)
        value = norm(phi(candidate)-candidate)/max(norm(candidate),eps);
    end

    function [value,isterminal,direction] = local_rate_event(~,candidate)
        value = rateLimitHz-max(abs(candidate));
        isterminal = 1;
        direction = -1;
    end
end
