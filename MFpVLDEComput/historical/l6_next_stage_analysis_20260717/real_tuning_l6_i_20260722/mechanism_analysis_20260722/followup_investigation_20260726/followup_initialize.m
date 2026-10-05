function paths = followup_initialize()
% Initialize the isolated follow-up investigation without changing old results.

root = fileparts(mfilename('fullpath'));
parent = fileparts(root);
addpath(parent,'-begin');
setupOverride = getenv('FOLLOWUP_SETUP_FILE');
if isempty(setupOverride)
    paths = mechanism_initialize();
    paths.ContinuationFile = fullfile(paths.OutputRoot,'branch_continuation.mat');
    paths.FollowupOutput = fullfile(root,'results');
else
    paths = struct();
    paths.SetupFile = setupOverride;
    paths.ContinuationFile = getenv('FOLLOWUP_CONTINUATION_FILE');
    paths.FollowupOutput = getenv('FOLLOWUP_OUTPUT_ROOT');
    if isempty(paths.ContinuationFile) || ~isfile(paths.ContinuationFile)
        error('Followup:ContinuationFile','FOLLOWUP_CONTINUATION_FILE is invalid.');
    end
    if isempty(paths.FollowupOutput)
        error('Followup:OutputRoot','FOLLOWUP_OUTPUT_ROOT is required.');
    end
end
addpath(root,'-begin');
paths.FollowupRoot = root;
paths.FollowupFigures = fullfile(paths.FollowupOutput,'figures');
if ~exist(paths.FollowupOutput,'dir'); mkdir(paths.FollowupOutput); end
if ~exist(paths.FollowupFigures,'dir'); mkdir(paths.FollowupFigures); end
end
