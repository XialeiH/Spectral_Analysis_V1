function repro_addpath(varargin)
% Keep study-specific helper versions isolated, preserving original path order.
root = fileparts(fileparts(mfilename('fullpath')));
workspace = fullfile(root,'Outputs','workspace');
for k=1:numel(varargin)
    value = char(varargin{k});
    if startsWith(value,[workspace filesep])
        relative = value(numel(workspace)+2:end);
        code = fullfile(root,'Utils',relative);
        if isfolder(code), varargin{k}=code; end
    end
end
addpath(varargin{:});
end
