function cfg = GA_h96_parameter_gain_config_from_env()
% Environment-driven multiplicative P perturbation for finite differences.
% GEOM_P_DIRECTION can be a 1-based index or one of GA_h96_parameter_spec names.
% GEOM_P_LOG_GAIN is the log gain, e.g. +0.01 or -0.01.
spec = GA_h96_parameter_spec();
cfg = struct();
cfg.Spec = spec;
cfg.Active = false;
cfg.Direction = [];
cfg.DirectionName = '';
cfg.LogGain = 0;
cfg.Gain = 1;
cfg.MatrixGains = ones(1, spec.NumParameters - 1);
cfg.L6KernelGain = 1;
cfg.Tag = 'baseline';

dirText = strtrim(getenv('GEOM_P_DIRECTION'));
if isempty(dirText)
    return;
end
logGain = str2double(getenv('GEOM_P_LOG_GAIN'));
if isnan(logGain)
    error('GEOM_P_LOG_GAIN must be numeric when GEOM_P_DIRECTION is set.');
end

dirIdx = str2double(dirText);
if isnan(dirIdx)
    dirIdx = find(strcmp(spec.Names, dirText), 1);
end
if isempty(dirIdx) || dirIdx < 1 || dirIdx > spec.NumParameters
    error('Unknown GEOM_P_DIRECTION %s.', dirText);
end

gain = exp(logGain);
cfg.Active = true;
cfg.Direction = dirIdx;
cfg.DirectionName = spec.Names{dirIdx};
cfg.LogGain = logGain;
cfg.Gain = gain;
if dirIdx <= numel(cfg.MatrixGains)
    cfg.MatrixGains(dirIdx) = gain;
else
    cfg.L6KernelGain = gain;
end

if logGain >= 0
    signText = 'plus';
else
    signText = 'minus';
end
cfg.Tag = sprintf('P%02d_%s_eps%s', dirIdx, cfg.DirectionName, signText);
cfg.Tag = regexprep(cfg.Tag, '[^A-Za-z0-9_]+', '_');
end
