
% Find the fixed point for ODE

function [value,isterminal,direction] = stopWhenFixed(t,f,ParaODE,tol_fixed)
    Phi_f = Phi_from_2D_sigmoid(f, ParaODE);
    rel   = norm(Phi_f - f, 2) / max(1e-12, norm(f,2));
    value = rel - tol_fixed;
    isterminal = 1;
    direction = -1;
end
