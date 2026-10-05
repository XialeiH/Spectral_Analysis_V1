function y = rescaledSigmoid(x, xi, s)
    % Validate input
    if s < 0 || s > 1
        error('Shape parameter must be between 0 and 1');
    end
    if xi < 0 || xi > 1
        error('Intersection point must be between 0 and 1');
    end

    % Adjust k based on s
    k = 10 * (0.01 + s * 9);  % Scale factor for steepness

    % Calculate the value of L to ensure that the curve passes through (xi, xi)
    L = (1 / xi) - 1;

    % Apply the general sigmoid function
    y = L ./ (1 + exp(-k * (x - xi)));

    % Normalize to ensure the curve passes through (0,0) and (1,1)
    min_y = L / (1 + exp(k * xi));
    max_y = L / (1 + exp(-k * (1 - xi)));
    y = (y - min_y) / (max_y - min_y);
end