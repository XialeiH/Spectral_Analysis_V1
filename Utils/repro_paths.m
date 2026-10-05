function location = repro_paths(relative)
% Resolve a retained input location inside the writable reproduction workspace.
root = fileparts(fileparts(mfilename('fullpath')));
location = fullfile(root,'Outputs','workspace',relative);
end
