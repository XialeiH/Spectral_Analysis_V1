function cachePath = build_h96_model_cache(runtimeDir, outputDir)
% Extract generated predictor constants once so iteration never reparses them.
if nargin < 2 || isempty(outputDir)
    outputDir = fileparts(mfilename('fullpath'));
end

populations = {'S','C','I'};
parameterNames = {'mu','sd','W1','b1','W2','b2','W3','b3','W4','b4'};
models = struct();

for p = 1:numel(populations)
    population = populations{p};
    sourcePath = fullfile(runtimeDir, sprintf('predict_pref6D_%s.m', population));
    assert(isfile(sourcePath), 'Missing predictor: %s', sourcePath);
    sourceText = fileread(sourcePath);
    model = struct();

    for k = 1:numel(parameterNames)
        name = parameterNames{k};
        token = regexp(sourceText, ['(?s)(?m)^' name '\s*=\s*(\[.*?\]);'], 'tokens', 'once');
        assert(~isempty(token), 'Could not extract %s from %s.', name, sourcePath);
        value = eval(token{1}); %#ok<EVLDIR>
        model.(name) = value;
    end
    models.(population) = model;
end

cachePath = fullfile(outputDir, 'h96_models_double.mat');
save(cachePath, 'models', '-v7');
end
