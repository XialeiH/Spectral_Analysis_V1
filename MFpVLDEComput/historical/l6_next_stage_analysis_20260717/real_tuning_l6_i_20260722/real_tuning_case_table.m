function cases = real_tuning_case_table(smoke)
% Define beta sweeps for whole-pathway real L6 and inhibition tuning.

if nargin < 1
    smoke = false;
end
if smoke
    betaGrid = [-0.50 -0.30 -0.20 -0.10 0 0.10 0.20 0.30 0.50]';
else
    betaGrid = (-0.50:0.025:0.50)';
end
representativeBeta = [-0.40 -0.20 0 0.20 0.40];

pathway = [repmat("L6",numel(betaGrid),1); ...
    repmat("Inhibition",numel(betaGrid),1)];
beta = [betaGrid;betaGrid];
isRepresentative = false(size(beta));
if ~smoke
    for index = 1:numel(representativeBeta)
        isRepresentative = isRepresentative | ...
            abs(beta-representativeBeta(index)) < 1e-12;
    end
end
taskId = (1:numel(beta))';
cases = table(taskId,pathway,beta,1+beta,isRepresentative, ...
    'VariableNames',{'taskId','pathway','beta','pathwayGain','isRepresentative'});
end
