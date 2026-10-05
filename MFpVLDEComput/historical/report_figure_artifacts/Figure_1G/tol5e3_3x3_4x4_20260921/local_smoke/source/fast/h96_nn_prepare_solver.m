function [fast, timing] = h96_nn_prepare_solver(ctx, contrast, orientation, precision, modelCachePath)
% Convenience preparation for a single h96 condition.
if nargin<4 || isempty(precision)
    precision = 'double';
end
if nargin<5
    modelCachePath = [];
end
timer = tic;
base = h96_nn_prepare_base(ctx, precision, modelCachePath);
timing.BaseSeconds = toc(timer);
timer = tic;
fast = h96_nn_set_condition(base, contrast, orientation);
timing.ConditionSeconds = toc(timer);
timing.TotalSeconds = timing.BaseSeconds+timing.ConditionSeconds;
end
