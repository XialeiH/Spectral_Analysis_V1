function sourceRoot=prepare_timing_source()
% Assemble the retained 3x3 runtime; this function submits no jobs.
root=fileparts(fileparts(mfilename('fullpath')));
sourceRoot=fullfile(root,'Outputs','timing_runtime');
response=fullfile(root,'Data','timing_runtime','bundles','cg_response_c100_a0.mat');
assert(isfile(response), ...
    'First run python3 Utils/restore_large_data.py from the repository root.');
if ~isfolder(sourceRoot), mkdir(sourceRoot); end
copyfile(fullfile(root,'Utils','timing_runtime','*'),sourceRoot,'f');
copyfile(fullfile(root,'Data','timing_runtime','*'),sourceRoot,'f');
end
