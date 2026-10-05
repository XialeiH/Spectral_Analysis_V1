function S = GA_structural_susceptibility(metricStack, paramVals)
% Finite-difference S_ija = partial g_ij / partial p_a.
% metricStack is nTask-by-nTask-by-nParamSample for one parameter coordinate.
paramVals = paramVals(:)';
nSample = numel(paramVals);
S = zeros(size(metricStack));

for k = 1:nSample
    if k == 1
        h = paramVals(2) - paramVals(1);
        S(:,:,k) = (metricStack(:,:,2) - metricStack(:,:,1)) ./ h;
    elseif k == nSample
        h = paramVals(end) - paramVals(end-1);
        S(:,:,k) = (metricStack(:,:,end) - metricStack(:,:,end-1)) ./ h;
    else
        h = paramVals(k+1) - paramVals(k-1);
        S(:,:,k) = (metricStack(:,:,k+1) - metricStack(:,:,k-1)) ./ h;
    end
end
end
