function row = mechanism_domain_stats(state,context,gain6,gainI,pathway,beta,stateLabel)
% Summarize trained-domain occupancy and extension effects at one state.

inputs = mechanism_state_inputs(state,context,gain6,gainI);
types = {'S','C','I'};
edgeE = 47500;
edgeI = 35315.63671875;
blendE = 5000;
blendI = 5000;

count = 0;
outside = 0;
blend = 0;
fullExtension = 0;
rawExtendedSq = 0;
extendedSq = 0;
responseAbove190 = 0;
gradientNormExtended = 0;
gradientNormRaw = 0;
maximumE = -Inf;
maximumI = -Inf;
for index = 1:numel(types)
    name = types{index};
    l4E = inputs.(name).L4E;
    l4I = inputs.(name).L4I;
    exE = max(0,l4E-edgeE);
    exI = max(0,l4I-edgeI);
    distance = sqrt((exE/blendE).^2+(exI/blendI).^2);
    outsideMask = distance>0;
    blendMask = distance>0 & distance<1;
    fullMask = distance>=1;
    [extended,gEE,gEI,gE6] = mechanism_aggregate_response( ...
        name,l4E,l4I,inputs.L6,context,'extended');
    [raw,gRE,gRI,gR6] = mechanism_aggregate_response( ...
        name,l4E,l4I,inputs.L6,context,'raw');
    count = count+numel(l4E);
    outside = outside+nnz(outsideMask);
    blend = blend+nnz(blendMask);
    fullExtension = fullExtension+nnz(fullMask);
    rawExtendedSq = rawExtendedSq+sum((extended-raw).^2);
    extendedSq = extendedSq+sum(extended.^2);
    responseAbove190 = responseAbove190+nnz(extended>190);
    gradientNormExtended = gradientNormExtended+sum(gEE.^2+gEI.^2+gE6.^2);
    gradientNormRaw = gradientNormRaw+sum(gRE.^2+gRI.^2+gR6.^2);
    maximumE = max(maximumE,max(l4E));
    maximumI = max(maximumI,max(l4I));
end

row = table(string(pathway),beta,string(stateLabel),gain6,gainI, ...
    outside/count,blend/count,fullExtension/count,responseAbove190/count, ...
    sqrt(rawExtendedSq)/max(sqrt(extendedSq),eps),maximumE,maximumI, ...
    min(inputs.L6),max(inputs.L6),sqrt(gradientNormExtended), ...
    sqrt(gradientNormRaw), ...
    'VariableNames',{'pathway','beta','stateLabel','gain6','gainI', ...
    'fractionOutsideTrainingBox','fractionInBlend','fractionFullExtension', ...
    'fractionResponseAbove190','rawVsExtendedRelativeResponseDifference', ...
    'maximumL4E','maximumL4I','minimumL6','maximumL6', ...
    'extendedAggregateGradientNorm','rawAggregateGradientNorm'});
end
