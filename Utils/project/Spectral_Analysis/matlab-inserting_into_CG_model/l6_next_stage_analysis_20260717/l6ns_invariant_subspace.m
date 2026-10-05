function subspace = l6ns_invariant_subspace(jacobian, count, selector, cfg)
% Partial Schur/Ritz invariant subspace obtained from converged eigenmodes.

modes = l6ns_eigenpairs(jacobian,count,selector,cfg);
q = orth(modes.Right);
projected = q'*(jacobian*q);
residual = norm(jacobian*q-q*projected,'fro')/max(norm(jacobian*q,'fro'),eps);

subspace = struct();
subspace.Basis = q;
subspace.RitzMatrix = projected;
subspace.RitzValues = eig(projected,'vector');
subspace.Residual = residual;
subspace.Modes = modes;
end
