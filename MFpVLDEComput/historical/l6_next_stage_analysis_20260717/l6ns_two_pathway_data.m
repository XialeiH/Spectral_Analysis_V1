function pathway = l6ns_two_pathway_data(data)
% Split the baseline Jacobian into rest, dynamic-L6, and I-source pathways.

n = data.PopulationSize;
iColumns = 2*n+(1:n);

pathway = struct();
pathway.J6 = sparse(data.B);
pathway.JI = sparse(data.Dimension,data.Dimension);
pathway.JI(:,iColumns) = data.J0(:,iColumns);
pathway.JRest = sparse(data.J0-pathway.J6-pathway.JI);
pathway.JBaseline = sparse(data.J0);
pathway.IColumns = iColumns;
pathway.PopulationSize = n;

pathway.ReconstructionError = norm( ...
    pathway.JRest+pathway.J6+pathway.JI-pathway.JBaseline,'fro') / ...
    max(norm(pathway.JBaseline,'fro'),eps);
pathway.J6ISourceFraction = norm(pathway.J6(:,iColumns),'fro') / ...
    max(norm(pathway.J6,'fro'),eps);
pathway.JINonISourceFraction = norm(pathway.JI(:,1:2*n),'fro') / ...
    max(norm(pathway.JI,'fro'),eps);
end
