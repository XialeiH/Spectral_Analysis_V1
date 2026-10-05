function output = h96_nn_all_population_responses_staticfirst(recurrentE,recurrentI,l6Axis,fast)
% Evaluate the h96 MLP with its static layer-1 contribution cached.
n = fast.n;
rows = fast.mixRows;
dynamicE = high_l4e_guard(reshape(recurrentE,n,3));
dynamicI = reshape(recurrentI,n,3);

x = zeros(fast.activePairCount,3,3,fast.nnPrecision);
x(:,1,:) = reshape((cast(l6Axis(rows),fast.nnPrecision)-fast.mu(:,3,:))./ ...
    fast.sd(:,3,:),fast.activePairCount,1,3);
x(:,2,:) = reshape((cast(dynamicE(rows,:),fast.nnPrecision)-reshape(fast.mu(:,4,:),1,3))./ ...
    reshape(fast.sd(:,4,:),1,3),fast.activePairCount,1,3);
x(:,3,:) = reshape((cast(dynamicI(rows,:),fast.nnPrecision)-reshape(fast.mu(:,5,:),1,3))./ ...
    reshape(fast.sd(:,5,:),1,3),fast.activePairCount,1,3);

x = tanh(fast.staticLayer1+pagemtimes(x,fast.dynamicW1T));
x = tanh(pagemtimes(x,fast.W2T)+fast.b2);
x = tanh(pagemtimes(x,fast.W3T)+fast.b3);
z = pagemtimes(x,fast.W4T)+fast.b4;
y = max(z,0)+log1p(exp(-abs(z)));
output = double(fast.mixAggregate*reshape(y,fast.activePairCount,3));
end

function x = high_l4e_guard(x)
edge = 47500;
full = 55000;
u = min(max((x-edge)./(full-edge),0),1);
g = u.^3.*(10-15.*u+6.*u.^2);
x = (1-g).*x+g.*(edge+0.02.*(x-edge));
end
