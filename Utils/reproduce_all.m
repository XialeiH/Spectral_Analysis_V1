function results = reproduce_all(includeHistorical)
% Run every current main/supplementary Live Script in a clean base workspace.
% Call from a dedicated MATLAB session: each entry intentionally clears variables.
if nargin < 1, includeHistorical=false; end
root=fileparts(fileparts(mfilename('fullpath')));
previous=pwd;
cleanup=onCleanup(@()cd(previous));
cd(root);
entries=jsondecode(fileread(fullfile(root,'panels.json')));
results=struct('id',{},'passed',{},'error',{});
for k=1:numel(entries)
    entry=entries(k);
    if strcmp(entry.group,'Additional') && ~includeHistorical, continue; end
    restoredefaultpath;
    addpath(fullfile(root,'Utils'));
    set(groot,'defaultFigureVisible','off');
    fprintf('Reproducing %s\n',entry.id);
    try
        evalin('base',sprintf('run(''%s'');',entry.live_script));
        result=struct('id',entry.id,'passed',true,'error','');
    catch exception
        result=struct('id',entry.id,'passed',false, ...
            'error',getReport(exception,'extended','hyperlinks','off'));
        fprintf(2,'%s\n',result.error);
    end
    evalin('base','clear cleanup');
    close all force;
    results(end+1)=result; %#ok<AGROW>
end
fid=fopen(fullfile(root,'Outputs','reproduction_results.json'),'w');
fprintf(fid,'%s',jsonencode(results)); fclose(fid);
assert(all([results.passed]),'One or more figures failed; see Outputs/reproduction_results.json.');
end
