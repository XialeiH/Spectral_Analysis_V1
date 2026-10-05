function output = h96_nn_predict_pages_codegen(dynamicL6, dynamicE, dynamicI, ...
    staticPages, mu, sd, W1T, b1, W2T, b2, W3T, b3, W4T, b4, pixelWeights)
%#codegen
% Fixed-size kernel compiled for the three h96 population predictors.
n = 1600;
x = staticPages;
output = zeros(n, 3, 'like', staticPages);

for population = 1:3
    x(:,3,population) = (dynamicL6 - mu(1,3,population)) ./ sd(1,3,population);
    x(:,4,population) = (dynamicE(:,population) - mu(1,4,population)) ./ sd(1,4,population);
    x(:,5,population) = (dynamicI(:,population) - mu(1,5,population)) ./ sd(1,5,population);

    activation = tanh(x(:,:,population) * W1T(:,:,population) + b1(:,:,population));
    activation = tanh(activation * W2T(:,:,population) + b2(:,:,population));
    activation = tanh(activation * W3T(:,:,population) + b3(:,:,population));
    z = activation * W4T(:,:,population) + b4(:,:,population);
    y = max(z, 0) + log1p(exp(-abs(z)));
    y = reshape(y, n, 5);
    output(:,population) = sum(pixelWeights .* y, 2);
end
end
