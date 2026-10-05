% Run MATLAB Code Analyzer on every file in this analysis suite.

codeRoot=repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');
repro_addpath(codeRoot);
files=dir(fullfile(codeRoot,'*.m'));
issueCount=0;
for fileIndex=1:numel(files)
    fileName=fullfile(files(fileIndex).folder,files(fileIndex).name);
    issues=checkcode(fileName,'-id');
    fprintf('%s: %d Code Analyzer issue(s).\n',files(fileIndex).name,numel(issues));
    for issueIndex=1:numel(issues)
        fprintf('  line %d col %d id %s: %s\n',issues(issueIndex).line, ...
            issues(issueIndex).column(1),issues(issueIndex).id,issues(issueIndex).message);
    end
    issueCount=issueCount+numel(issues);
end
fprintf('L6NS_CHECKCODE_TOTAL=%d\n',issueCount);
