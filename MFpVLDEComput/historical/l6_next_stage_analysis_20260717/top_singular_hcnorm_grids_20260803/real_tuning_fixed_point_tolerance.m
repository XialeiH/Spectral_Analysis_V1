function result = real_tuning_fixed_point_tolerance(phi, initialState, relaxationP, tolerance)
% Compute a moved fixed point to the requested relative residual tolerance.

state = initialState(:);
maximumIterations = 2500;
rateLimitHz = 500;
termination = "maximum_iterations";
response = phi(state);
residual = Inf;

for iteration = 1:maximumIterations
    next = (1-relaxationP)*state + relaxationP*response;
    nextResponse = phi(next);
    residual = norm(nextResponse-next)/max(norm(next),eps);
    state = next;
    if any(~isfinite(state))
        termination = "nonfinite";
        break
    end
    if max(abs(state)) >= rateLimitHz
        termination = "rate_limit";
        break
    end
    if residual <= tolerance
        termination = "converged";
        break
    end
    response = nextResponse;
end

result = struct('State',state,'Residual',residual,'Iterations',iteration, ...
    'Converged',isfinite(residual) && residual<=tolerance, ...
    'Termination',termination,'MinimumRateHz',min(state), ...
    'MaximumRateHz',max(state));
end
