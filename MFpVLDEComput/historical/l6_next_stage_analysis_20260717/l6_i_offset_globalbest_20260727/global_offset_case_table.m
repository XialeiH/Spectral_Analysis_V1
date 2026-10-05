function cases = global_offset_case_table(smoke)
% L6 gain grid for globally audited inhibitory-offset optimization.

if nargin<1; smoke=false; end
if smoke
    beta6 = [0 0.04 0.25]';
else
    beta6 = (0:0.002:0.25)';
end
taskId = (1:numel(beta6))';
cases = table(taskId,beta6);
end
