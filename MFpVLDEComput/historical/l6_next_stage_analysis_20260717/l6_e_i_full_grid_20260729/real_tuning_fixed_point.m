function result = real_tuning_fixed_point(phi, initialState, relaxationP)
% Compute the moved fixed point using the model's regular relaxed iteration.

state = initialState(:);
maximumIterations = 2500;
tolerance = 1e-10;
rateLimitHz = 500;
residualHistory = nan(maximumIterations,1);
stepHistory = nan(maximumIterations,1);
termination = "maximum_iterations";

response = phi(state);
for iteration = 1:maximumIterations
    next = (1-relaxationP)*state + relaxationP*response;
    nextResponse = phi(next);
    residualHistory(iteration) = norm(nextResponse-next)/max(norm(next),eps);
    stepHistory(iteration) = norm(next-state)/max(norm(next),eps);
    state = next;
    if any(~isfinite(state))
        termination = "nonfinite";
        break
    end
    if max(abs(state)) >= rateLimitHz
        termination = "rate_limit";
        break
    end
    if residualHistory(iteration) <= tolerance
        termination = "converged";
        break
    end
    response = nextResponse;
end

residualHistory = residualHistory(1:iteration);
stepHistory = stepHistory(1:iteration);
result = struct('State',state,'Residual',residualHistory(end), ...
    'Iterations',iteration,'Converged',residualHistory(end)<=1e-9, ...
    'Termination',termination,'RateLimitHz',rateLimitHz, ...
    'MinimumRateHz',min(state),'MaximumRateHz',max(state), ...
    'ResidualHistory',residualHistory,'StepHistory',stepHistory);
end
