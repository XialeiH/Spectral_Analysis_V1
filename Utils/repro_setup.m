function state = repro_setup(id)
% Copy immutable inputs once; every generated file stays under Outputs.
root = fileparts(fileparts(mfilename('fullpath')));
dataRoot = fullfile(root,'Data');
workspace = fullfile(root,'Outputs','workspace');
revision = strtrim(fileread(fullfile(dataRoot,'REVISION')));
marker = fullfile(workspace,'REVISION');
if ~isfile(marker) || ~strcmp(strtrim(fileread(marker)),revision)
    if ~isfolder(workspace), mkdir(workspace); end
    groups = {'project','benchmarks'};
    for k=1:numel(groups)
        source = fullfile(dataRoot,groups{k});
        if isfolder(source), copyfile(source,fullfile(workspace,groups{k}),'f'); end
    end
    copyfile(fullfile(dataRoot,'REVISION'),marker,'f');
end
state = struct('ID',id,'Started',now,'Root',root,'Workspace',workspace);
state.Before=dir(fullfile(workspace,'**','*'));
state.Before=state.Before(~[state.Before.isdir]);
end
