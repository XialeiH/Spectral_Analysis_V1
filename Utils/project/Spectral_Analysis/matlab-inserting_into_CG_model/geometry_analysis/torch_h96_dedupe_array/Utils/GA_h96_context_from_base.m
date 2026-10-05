function ctx = GA_h96_context_from_base()
% Collect the standard h96 workspace variables needed by the response map.
names = {'PixLGNCtgr','L6Kernel', ...
    'C_SS_meanU','C_CS_meanU','C_IS_mean', ...
    'C_SC_meanU','C_CC_meanU','C_IC_mean', ...
    'C_SI_mean','C_CI_mean','C_II_mean', ...
    'L4SEp','L4SIp','L4CEp','L4CIp','L4IEp','L4IIp', ...
    'L4EmeshXAll','L4ImeshYAll','LDEFrfuncAll', ...
    'N_HCOutY','NPixX','NPixY','Isaturation','EKpUse','IKpUse', ...
    'L6pars','L6parId','p','LDEsigmoid','LDEEquv','IniTest'};

for k = 1:numel(names)
    name = names{k};
    if evalin('base', sprintf('exist(''%s'', ''var'')', name))
        ctx.(name) = evalin('base', name);
    end
end

required = {'PixLGNCtgr','L6Kernel','C_SS_meanU','C_CS_meanU','C_IS_mean', ...
    'C_SC_meanU','C_CC_meanU','C_IC_mean','C_SI_mean','C_CI_mean','C_II_mean', ...
    'L4SEp','L4SIp','L4CEp','L4CIp','L4IEp','L4IIp', ...
    'N_HCOutY','NPixX','NPixY','Isaturation','EKpUse','IKpUse','L6pars','L6parId'};
missing = {};
for k = 1:numel(required)
    if ~isfield(ctx, required{k})
        missing{end+1} = required{k}; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('GA_h96_context_from_base:MissingVariables', ...
        'Missing h96 setup variables in base workspace: %s', strjoin(missing, ', '));
end

if ~isfield(ctx, 'p')
    ctx.p = 0.33;
end
if ~isfield(ctx, 'L4EmeshXAll')
    ctx.L4EmeshXAll = {[]};
end
if ~isfield(ctx, 'L4ImeshYAll')
    ctx.L4ImeshYAll = {[]};
end
if ~isfield(ctx, 'LDEFrfuncAll')
    ctx.LDEFrfuncAll = {struct('S',{{}},'C',{{}},'I',{{}})};
end
end
