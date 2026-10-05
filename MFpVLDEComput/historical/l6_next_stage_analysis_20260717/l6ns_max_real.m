function [value, lambda, rightVector] = l6ns_max_real(jacobian, cfg)
% Largest real part of the Jacobian spectrum.

modes = l6ns_eigenpairs(jacobian, min(6, size(jacobian, 1) - 2), 'largestreal', cfg);
[value, index] = max(real(modes.Lambda));
lambda = modes.Lambda(index);
rightVector = modes.Right(:, index);
end
