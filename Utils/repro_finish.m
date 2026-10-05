function outputs = repro_finish(state)
% Gather newly exported panels without publishing the input workspace.
destination = fullfile(state.Root,'Outputs',state.ID);
if ~isfolder(destination), mkdir(destination); end
outputs = {};
for extension = {'.pdf','.fig','.csv','.tsv'}
    files = dir(fullfile(state.Workspace,'**',['*' extension{1}]));
    for k=1:numel(files)
        source = fullfile(files(k).folder,files(k).name);
        if isfield(state,'Before')
            prior=find(strcmp({state.Before.folder},files(k).folder) & ...
                strcmp({state.Before.name},files(k).name),1);
            changed=isempty(prior) || files(k).datenum~=state.Before(prior).datenum || ...
                files(k).bytes~=state.Before(prior).bytes;
        else
            changed=files(k).datenum>=state.Started;
        end
        if changed
            name = files(k).name;
            if startsWith(state.ID,'Main_')
                [~,stem,ext]=fileparts(name);
                if strcmp(state.ID,'Main_04_B') && strcmp(stem,'Figure 4D')
                    continue
                end
                mapping=struct('Main_01_B','Figure 1B_only_iterations', ...
                    'Main_01_C','Figure 1B','Main_01_D','Figure 1C', ...
                    'Main_01_E','Figure 1D','Main_01_FG','Figure 1E', ...
                    'Main_02_B','Figure 2B','Main_02_C','Figure 2D', ...
                    'Main_02_D','Figure 2E','Main_04_A','Figure 3C', ...
                    'Main_04_B','Figure 3D','Main_04_C','Figure 3E.2', ...
                    'Main_05_A','Figure 4A','Main_05_B','Figure 4B', ...
                    'Main_05_C','Figure 4C','Main_05_D','Figure 4E', ...
                    'Main_05_E','Figure 4F','Main_06_A','Figure 5A', ...
                    'Main_06_B','Figure 5B','Main_06_C','Figure 5C', ...
                    'Main_06_D','Figure 5D','Main_06_E','Figure 5E.2', ...
                    'Main_06_F','Figure 5E.1');
                if isfield(mapping,state.ID) && startsWith(stem,mapping.(state.ID))
                    prefix=mapping.(state.ID);
                    name=[state.ID stem(numel(prefix)+1:end) ext];
                elseif strcmp(state.ID,'Main_03_AB') && strcmp(stem,'Figure 2C.1')
                    name=['Main_03_A' ext];
                elseif strcmp(state.ID,'Main_03_AB') && strcmp(stem,'Figure 2C.2')
                    name=['Main_03_B' ext];
                elseif strcmp(state.ID,'Main_04_B') && strcmp(stem,'Figure_4D_summary')
                    name=['Main_04_B_summary' ext];
                end
            end
            target = fullfile(destination,name);
            copyfile(source,target,'f');
            outputs{end+1,1} = target; %#ok<AGROW>
        end
    end
end
assert(~isempty(outputs),'No new figure exports were produced.');
disp(outputs);
end
