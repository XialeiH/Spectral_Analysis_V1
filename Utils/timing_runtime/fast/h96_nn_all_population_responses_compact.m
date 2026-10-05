function output = h96_nn_all_population_responses_compact(recurrentE, recurrentI, l6Axis, fast)
% Evaluate only pixel/LGN pairs having nonzero mixture weight.
n = fast.n;
rows = fast.mixRows;
dynamicE = h96_high_l4e_guard(reshape(recurrentE, n, 3));
dynamicI = reshape(recurrentI, n, 3);

x = fast.staticCompactPages;
x(:,3,:) = reshape((cast(l6Axis(rows), fast.nnPrecision)-fast.mu(:,3,:))./ ...
    fast.sd(:,3,:), fast.activePairCount, 1, 3);
x(:,4,:) = reshape((cast(dynamicE(rows,:), fast.nnPrecision)-reshape(fast.mu(:,4,:), 1, 3))./ ...
    reshape(fast.sd(:,4,:), 1, 3), fast.activePairCount, 1, 3);
x(:,5,:) = reshape((cast(dynamicI(rows,:), fast.nnPrecision)-reshape(fast.mu(:,5,:), 1, 3))./ ...
    reshape(fast.sd(:,5,:), 1, 3), fast.activePairCount, 1, 3);

x = tanh(pagemtimes(x, fast.W1T)+fast.b1);
x = tanh(pagemtimes(x, fast.W2T)+fast.b2);
x = tanh(pagemtimes(x, fast.W3T)+fast.b3);
z = pagemtimes(x, fast.W4T)+fast.b4;
y = max(z, 0)+log1p(exp(-abs(z)));
y = reshape(y, fast.activePairCount, 3);
output = double(fast.mixAggregate*y);
end

function xg = h96_high_l4e_guard(x)
edge = 47500;
full = 55000;
tailSlope = 0.02;
u = min(max((x-edge)./(full-edge), 0), 1);
g = u.^3.*(10-15.*u+6.*u.^2);
xg = (1-g).*x+g.*(edge+tailSlope.*(x-edge));
end
