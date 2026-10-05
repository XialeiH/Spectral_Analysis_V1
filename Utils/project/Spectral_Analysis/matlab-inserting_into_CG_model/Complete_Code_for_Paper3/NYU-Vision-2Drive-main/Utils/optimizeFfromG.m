function f = optimizeFfromG(g, aa)
    % Validate inputs
    if aa <= 0 || 180 / aa ~= round(180 / aa)
        error('aa must be a positive number that evenly divides 180');
    end
    
    % g is the given vector where g(x) = f(x) + f(x-45), x = 0:aa:180
    g = g(:); % Ensure g is a column vector for consistency
    N = length(g); % Number of points in g

    % Calculate the shift for x-45 based on the discretization step aa
    shift = 45 / aa;

    % Initial guess for f, could be based on g for simplicity
    f0 = max(0, g/2); % Start with a positive assumption

    % Define the symmetry constraints dynamically based on aa
    midPoint = N/2 + 0.5; % Adjust if N is odd
    Aeq = [];
    beq = [];
    for i = 1:floor(N/2) % Ensure symmetry
        row = zeros(1, N);
        row(i) = 1;
        row(N-i+1) = -1;
        Aeq = [Aeq; row];
        beq = [beq; 0];
    end

    % Bounds (f(x) >= 0)
    lb = zeros(N, 1);
    ub = []; % No upper bounds
    
    % Options for fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

    % Objective function: Minimize the difference between actual and modeled g(x)
    objective = @(f) sum((g - modeledG(f, shift, N)).^2);

    % Solve using fmincon
    [f, ~] = fmincon(objective, f0, [], [], Aeq, beq, lb, ub, @nonlcon, options);

    % Modeled g(x) based on current f, considering the shift for f(x-45)
    function g_mod = modeledG(f, shift, N)
        f_extended = [f(1:end-1); f]; % Extend to handle periodic boundary
        % Calculate indices for f(x) and f(x-45)
        indices = 1:N;
        shiftIndices = mod(indices - 1 - shift, N) + 1;
        g_mod = (f + f_extended(shiftIndices))/2;
    end

    % Nonlinear constraints: single peak, min at 0, max at 90
    function [c, ceq] = nonlcon(f)
        % Approximate derivative to check for single peak
        df = diff(f);
        % Ensure non-increasing on both sides of the peak at 90 degrees
        peakIndex = round(90 / aa) + 1;
        c = [df(1:peakIndex-1); -df(peakIndex:end)];
        % No equality constraints
        ceq = [];
    end
end
