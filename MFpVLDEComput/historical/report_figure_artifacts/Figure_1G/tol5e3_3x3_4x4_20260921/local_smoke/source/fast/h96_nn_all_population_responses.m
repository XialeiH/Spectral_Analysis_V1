function output = h96_nn_all_population_responses(recurrentE, recurrentI, l6Axis, fast)
% Evaluate S/C/I h96 predictors together using page-wise matrix products.
n = fast.n;
dynamicL6 = repmat(l6Axis, 5, 1);
dynamicE = repmat(h96_high_l4e_guard(reshape(recurrentE, n, 3)), 5, 1);
dynamicI = repmat(reshape(recurrentI, n, 3), 5, 1);

x = fast.staticPages;
x(:,3,:) = reshape((cast(dynamicL6, fast.nnPrecision)-fast.mu(:,3,:))./ ...
    fast.sd(:,3,:), n*5, 1, 3);
x(:,4,:) = reshape((cast(dynamicE, fast.nnPrecision)-reshape(fast.mu(:,4,:), 1, 3))./ ...
    reshape(fast.sd(:,4,:), 1, 3), n*5, 1, 3);
x(:,5,:) = reshape((cast(dynamicI, fast.nnPrecision)-reshape(fast.mu(:,5,:), 1, 3))./ ...
    reshape(fast.sd(:,5,:), 1, 3), n*5, 1, 3);

x = tanh(pagemtimes(x, fast.W1T)+fast.b1);
x = tanh(pagemtimes(x, fast.W2T)+fast.b2);
x = tanh(pagemtimes(x, fast.W3T)+fast.b3);
z = pagemtimes(x, fast.W4T)+fast.b4;
y = max(z, 0)+log1p(exp(-abs(z)));
y = reshape(y, n, 5, 3);
output = double(reshape(sum(fast.pixelWeights.*y, 2), n, 3));
end

function xg = h96_high_l4e_guard(x)
edge = 47500;
full = 55000;
tailSlope = 0.02;
u = min(max((x-edge)./(full-edge), 0), 1);
g = u.^3.*(10-15.*u+6.*u.^2);
xg = (1-g).*x+g.*(edge+tailSlope.*(x-edge));
end
