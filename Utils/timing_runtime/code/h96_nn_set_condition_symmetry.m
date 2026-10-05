function fast = h96_nn_set_condition_symmetry(base,contrast,orientation)
% Cache layer-1 constants and the repeated 2-by-2-HC spatial classes.
fast = h96_nn_set_condition_staticfirst(base,contrast,orientation);
[row,column] = ind2sub([fast.fieldRows,fast.fieldCols],fast.mixRows);
periodRows = 20;
periodColumns = 20;
key = [mod(row-1,periodRows),mod(column-1,periodColumns),fast.mixCategories];
[~,fast.symRepresentative,fast.symGroup] = unique(key,'rows');
fast.symClassCount = numel(fast.symRepresentative);
end
