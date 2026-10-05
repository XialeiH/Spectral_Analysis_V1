function run_leading_mode_component_audit()
% Attribute leading eigenvalues to E, L6, and inhibitory-source components.

paths = mechanism_initialize();
loaded = load(paths.SetupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
summary = readtable(fullfile(paths.ResultsRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');
pathways = ["L6";"L6";"L6";"Inhibition";"Inhibition"];
betas = [.15;.175;.20;-.05;-.075];
rows = cell(2*numel(betas),1);

for index = 1:numel(betas)
    pathway = pathways(index);
    beta = betas(index);
    moved = local_saved_state(paths.ResultsRoot,summary,pathway,beta);
    gain6 = 1;
    gainI = 1;
    if pathway=="L6"; gain6 = 1+beta; else; gainI = 1+beta; end
    rows{2*index-1} = local_one_state( ...
        baseline,context,gain6,gainI,pathway,beta,"baseline state");
    rows{2*index} = local_one_state( ...
        moved,context,gain6,gainI,pathway,beta,"moved equilibrium");
end

componentTable = vertcat(rows{:});
writetable(componentTable,fullfile(paths.OutputRoot, ...
    'leading_mode_component_attribution.tsv'),'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'leading_mode_component_attribution.mat'), ...
    'componentTable','-v7.3');
end

function row = local_one_state(state,context,gain6,gainI,pathway,beta,stateLabel)
[J,components,details] = mechanism_jacobian_components(state,context,gain6,gainI);
jReference = real_tuning_true_jacobian(state,context,gain6,gainI);
reconstructionError = norm(J-jReference,'fro')/max(norm(jReference,'fro'),eps);
[lambda,right,left] = local_leading_pair(J);
denominator = left'*right;
cE = (left'*(components.Excitatory*right))/denominator;
c6 = (left'*(components.L6*right))/denominator;
cI = (left'*(components.Inhibition*right))/denominator;
row = table(string(pathway),beta,string(stateLabel),real(lambda),imag(lambda), ...
    real(cE),real(c6),real(cI),abs(lambda-(cE+c6+cI)), ...
    norm(components.Excitatory,'fro'),norm(components.L6,'fro'), ...
    norm(components.Inhibition,'fro'),reconstructionError, ...
    median(details.EAdjustment),median(details.IGradient), ...
    median(abs(details.L4EGradient)),median(abs(details.L4IGradient)), ...
    median(abs(details.L6Gradient)), ...
    'VariableNames',{'pathway','beta','stateLabel','leadingReal','leadingImag', ...
    'leadingContributionExcitatory','leadingContributionL6', ...
    'leadingContributionInhibition','contributionReconstructionError', ...
    'frobeniusExcitatory','frobeniusL6','frobeniusInhibition', ...
    'jacobianReconstructionError','medianEAdjustment','medianIGradient', ...
    'medianAbsResponseGradientL4E','medianAbsResponseGradientL4I', ...
    'medianAbsResponseGradientL6'});
fprintf(['%s beta=%+.3f %s: lambda=%.6f, contributions E/L6/I = ' ...
    '%+.6f / %+.6f / %+.6f\n'],pathway,beta,stateLabel,real(lambda), ...
    real(cE),real(c6),real(cI));
end

function [lambda,right,left] = local_leading_pair(J)
options = struct('tol',1e-10,'maxit',2200,'p',100,'isreal',true,'disp',0);
[rightVectors,rightValues] = eigs(J,8,'largestreal',options);
rightValues = diag(rightValues);
[~,rightIndex] = max(real(rightValues));
lambda = rightValues(rightIndex);
right = rightVectors(:,rightIndex);
[leftVectors,leftValues] = eigs(J',12,'largestreal',options);
leftValues = diag(leftValues);
[~,leftIndex] = min(abs(leftValues-conj(lambda)));
left = leftVectors(:,leftIndex);
end

function state = local_saved_state(resultsRoot,summary,pathway,beta)
mask = summary.pathway==pathway & abs(summary.beta-beta)<1e-10;
taskId = summary.taskId(mask);
match = dir(fullfile(resultsRoot,sprintf('%03d_*',taskId),'case_result.mat'));
loaded = load(fullfile(match.folder,match.name),'moved');
state = loaded.moved(:);
end
