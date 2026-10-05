function cases = fpp_offset_case_table(smoke)
% L6 increases used for the FPP constant-stability offset experiment.

if nargin<1; smoke=false; end
if smoke
    beta6 = [0 0.02 0.04]';
else
    beta6 = (0:0.002:0.04)';
end
taskId = (1:numel(beta6))';
cases = table(taskId,beta6);
end
