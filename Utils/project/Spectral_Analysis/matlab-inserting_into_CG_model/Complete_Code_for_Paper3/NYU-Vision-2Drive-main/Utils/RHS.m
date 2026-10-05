% ODE is defined here

function dfdt = RHS(t, f, ParaODE)
    Phi_f = Phi_from_2D_sigmoid(f, ParaODE);
    dfdt  = (-f + Phi_f) / ParaODE.tau;
end