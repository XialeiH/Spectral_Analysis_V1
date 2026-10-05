function audit_revised_jfull_cluster_data(originalRoot, revisedRoot)
% Verify that only the four intended J_full cluster fields changed.

if nargin < 1 || isempty(originalRoot)
    originalRoot = getenv('REVISED_CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 2 || isempty(revisedRoot)
    revisedRoot = getenv('REVISED_CLUSTER_OUTPUT_ROOT');
end
files = dir(fullfile(originalRoot, 'full_eigenspectrum_maps_angle*_contrast*.mat'));
clusterFields = {'TopClusterEnvelopeMaps', 'TopClusterEigenvalues', ...
    'TopClusterCount', 'TopClusterRank'};
if numel(files) ~= 16
    error('RevisedClusterAudit:FileCount', 'Expected 16 original condition files.');
end
for fileIndex = 1:numel(files)
    original = load(fullfile(originalRoot, files(fileIndex).name), 'conditionData');
    revised = load(fullfile(revisedRoot, files(fileIndex).name), 'conditionData');
    originalAnalysis = rmfield(original.conditionData.Jacobians.J_full, clusterFields);
    revisedAnalysis = rmfield(revised.conditionData.Jacobians.J_full, clusterFields);
    original.conditionData.Jacobians.J_full = originalAnalysis;
    revised.conditionData.Jacobians.J_full = revisedAnalysis;
    if ~isequaln(original.conditionData, revised.conditionData)
        error('RevisedClusterAudit:UnexpectedChange', ...
            'Unexpected non-cluster change in %s.', files(fileIndex).name);
    end
end
fprintf('PASS: all 16 files differ only in the four J_full TopCluster fields.\n');
end
