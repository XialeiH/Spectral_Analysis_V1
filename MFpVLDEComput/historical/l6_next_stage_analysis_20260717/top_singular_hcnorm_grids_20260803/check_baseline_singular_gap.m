function check_baseline_singular_gap()
% Verify that the baseline leading right singular vector is identifiable.

setupFile = getenv('TSHC_SETUP_FILE');
outputRoot = getenv('TSHC_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('TopSingularHC:Environment','Setup file and output root are required.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(setupFile,'setup');
J = loaded.setup.Pathway.JBaseline;
options = struct('tol',1e-9,'maxit',2500,'p',80,'disp',0);
timer = tic;
[leftVectors,singularValues,rightVectors] = svds(J,2,'largest',options);
elapsedSeconds = toc(timer);
singularValues = diag(singularValues);
[singularValues,order] = sort(real(singularValues),'descend');
leftVectors = leftVectors(:,order);
rightVectors = rightVectors(:,order);
[rightMode,leftMode] = local_real_singular_pair( ...
    rightVectors(:,1),leftVectors(:,1));
gapAbsolute = singularValues(1)-singularValues(2);
gapRelative = gapAbsolute/max(singularValues(1),eps);
rightResidual = norm(J*rightMode-singularValues(1)*leftMode)/ ...
    max(singularValues(1),eps);
leftResidual = norm(J'*leftMode-singularValues(1)*rightMode)/ ...
    max(singularValues(1),eps);

summary = table(singularValues(1),singularValues(2),gapAbsolute, ...
    gapRelative,rightResidual,leftResidual,elapsedSeconds, ...
    'VariableNames',{'sigma1','sigma2','gapAbsolute','gapRelative', ...
    'rightResidual','leftResidual','elapsedSeconds'});
writetable(summary,fullfile(outputRoot,'baseline_singular_gap.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'baseline_top_singular_mode.mat'), ...
    'rightMode','leftMode','singularValues','summary','-v7.3');
disp(summary);
end

function [rightVector,leftVector] = local_real_singular_pair(rightVector,leftVector)
[~,index] = max(abs(rightVector));
phase = exp(-1i*angle(rightVector(index)));
rightVector = real(rightVector*phase);
leftVector = real(leftVector*phase);
rightVector = rightVector/max(norm(rightVector),eps);
leftVector = leftVector/max(norm(leftVector),eps);
end
