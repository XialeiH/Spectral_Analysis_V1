function output = h96_nn_all_population_responses_cached(recurrentE,recurrentI,l6Axis,fast,cache)
% Reuse exact spatial equivalence classes discovered at the first epoch.
n = fast.n;
rows = fast.mixRows;
dynamicE = high_l4e_guard(reshape(recurrentE,n,3));
dynamicI = reshape(recurrentI,n,3);
l6Rows = cast(l6Axis(rows),fast.nnPrecision);
output = zeros(n,3);
for population = 1:3
    representative = cache{population}.Representative;
    group = cache{population}.Group;
    dynamic = zeros(fast.activePairCount,3,fast.nnPrecision);
    dynamic(:,1) = (l6Rows-fast.mu(1,3,population))./fast.sd(1,3,population);
    dynamic(:,2) = (cast(dynamicE(rows,population),fast.nnPrecision)- ...
        fast.mu(1,4,population))./fast.sd(1,4,population);
    dynamic(:,3) = (cast(dynamicI(rows,population),fast.nnPrecision)- ...
        fast.mu(1,5,population))./fast.sd(1,5,population);
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
