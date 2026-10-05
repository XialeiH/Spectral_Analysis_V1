function jacobian = l6ns_pathway_jacobian(pathway,gamma6,gammaI)
% Evaluate J_rest + gamma6 J_6 + gammaI J_I.

jacobian = pathway.JRest + gamma6*pathway.J6 + gammaI*pathway.JI;
end
