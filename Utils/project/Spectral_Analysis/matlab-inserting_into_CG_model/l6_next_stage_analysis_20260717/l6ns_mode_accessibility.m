function tableOut = l6ns_mode_accessibility(modes,operators,cfg)
% Controllability, observability, and lab-identifiability of leading modes.

count = numel(modes.Lambda);
rows = cell(count,1);
for k = 1:count
    r = modes.Right(:,k);
    l = modes.Left(:,k);
    controlCoupling = norm(l'*operators.B)/max(norm(l),eps);
    observationCoupling = norm(operators.C*r)/max(norm(r),eps);
    etaLab = controlCoupling*observationCoupling / ...
        max(abs(1-modes.Lambda(k)),eps);
    rows{k} = {k,real(modes.Lambda(k)),imag(modes.Lambda(k)), ...
        (real(modes.Lambda(k))-1)/cfg.TauMs,controlCoupling, ...
        observationCoupling,etaLab};
end
tableOut = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'mode','lambdaReal','lambdaImag','alphaPhysicalPerMs', ...
    'controlAccessibility','observability','labIdentifiability'});
end
