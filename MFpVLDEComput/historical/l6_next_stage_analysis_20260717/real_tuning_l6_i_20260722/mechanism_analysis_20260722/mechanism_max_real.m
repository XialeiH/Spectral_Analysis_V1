function value = mechanism_max_real(matrix)
% Largest real part of the map Jacobian spectrum.

options = struct('tol',1e-9,'maxit',1800,'p',80,'isreal',true,'disp',0);
try
    eigenvalues = eigs(matrix,8,'largestreal',options);
catch
    eigenvalues = eigs(matrix,8,'lr',options);
end
value = max(real(eigenvalues));
end
