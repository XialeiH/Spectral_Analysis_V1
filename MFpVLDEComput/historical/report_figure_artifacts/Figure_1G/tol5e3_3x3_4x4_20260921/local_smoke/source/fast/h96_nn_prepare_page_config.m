function page = h96_nn_prepare_page_config(fast, precision)
% Stack the S/C/I h96 predictors so all three evaluate page-wise together.
if nargin < 2
    precision = 'double';
end
assert(strcmp(precision, 'double') || strcmp(precision, 'single'), ...
    'Precision must be double or single.');

page = fast;
page.nnPrecision = precision;
names = {'S', 'C', 'I'};
nRows = fast.n * 5;
page.staticPages = zeros(nRows, 5, 3, precision);
page.mu = zeros(1, 5, 3, precision);
page.sd = zeros(1, 5, 3, precision);
page.W1T = zeros(5, 96, 3, precision);
page.b1 = zeros(1, 96, 3, precision);
page.W2T = zeros(96, 96, 3, precision);
page.b2 = zeros(1, 96, 3, precision);
page.W3T = zeros(96, 96, 3, precision);
page.b3 = zeros(1, 96, 3, precision);
page.W4T = zeros(96, 1, 3, precision);
page.b4 = zeros(1, 1, 3, precision);

for population = 1:3
    name = names{population};
    model = fast.models.(name);
    page.staticPages(:,:,population) = cast(fast.staticNormalized.(name), precision);
    page.mu(:,:,population) = cast(model.mu, precision);
    page.sd(:,:,population) = cast(model.sd, precision);
    page.W1T(:,:,population) = cast(model.W1.', precision);
    page.b1(:,:,population) = cast(model.b1, precision);
    page.W2T(:,:,population) = cast(model.W2.', precision);
    page.b2(:,:,population) = cast(model.b2, precision);
    page.W3T(:,:,population) = cast(model.W3.', precision);
    page.b3(:,:,population) = cast(model.b3, precision);
    page.W4T(:,:,population) = cast(model.W4.', precision);
    page.b4(:,:,population) = cast(model.b4, precision);
end

page.pixelWeights = cast(reshape(fast.PixLGNCtgr, fast.n, 5, 1), precision);
end
