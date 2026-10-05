function [output,uniqueCounts,cache] = h96_nn_all_population_responses_unique(recurrentE,recurrentI,l6Axis,fast)
% Evaluate repeated predictor tuples once, then expand to the full field.
n = fast.n;
rows = fast.mixRows;
categories = fast.mixCategories;
dynamicE = high_l4e_guard(reshape(recurrentE,n,3));
dynamicI = reshape(recurrentI,n,3);
l6Rows = cast(l6Axis(rows),fast.nnPrecision);
output = zeros(n,3);
uniqueCounts = zeros(1,3);
cache = cell(3,1);

for population = 1:3
    dynamic = zeros(fast.activePairCount,3,fast.nnPrecision);
    dynamic(:,1) = (l6Rows-fast.mu(1,3,population))./fast.sd(1,3,population);
    dynamic(:,2) = (cast(dynamicE(rows,population),fast.nnPrecision)- ...
        fast.mu(1,4,population))./fast.sd(1,4,population);
    dynamic(:,3) = (cast(dynamicI(rows,population),fast.nnPrecision)- ...
        fast.mu(1,5,population))./fast.sd(1,5,population);
    key = [cast(categories,fast.nnPrecision),dynamic];
    [~,representative,group] = unique(key,'rows');
    cache{population} = struct('Representative',representative,'Group',group);
    uniqueCounts(population) = numel(representative);
    h = tanh(fast.staticLayer1(representative,:,population)+ ...
        dynamic(representative,:)*fast.dynamicW1T(:,:,population));
    h = tanh(h*fast.W2T(:,:,population)+fast.b2(:,:,population));
    h = tanh(h*fast.W3T(:,:,population)+fast.b3(:,:,population));
    z = h*fast.W4T(:,:,population)+fast.b4(:,:,population);
    y = max(z,0)+log1p(exp(-abs(z)));
    output(:,population) = double(fast.mixAggregate*y(group));
end
end

function x = high_l4e_guard(x)
edge = 47500;
full = 55000;
u = min(max((x-edge)./(full-edge),0),1);
g = u.^3.*(10-15.*u+6.*u.^2);
x = (1-g).*x+g.*(edge+0.02.*(x-edge));
end
