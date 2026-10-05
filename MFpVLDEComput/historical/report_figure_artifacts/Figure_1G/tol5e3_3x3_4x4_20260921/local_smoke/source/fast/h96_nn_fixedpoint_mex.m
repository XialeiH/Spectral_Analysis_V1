function [state, info] = h96_nn_fixedpoint_mex(initialState, fast, iterations)
% Fixed-count h96 iteration using a generated fixed-size NN predictor MEX.
S = initialState.S;
C = initialState.C;
I = initialState.I;
n = fast.n;

timer = tic;
for epoch = 1:iterations
    E = 0.6923 .* S + 0.3077 .* C;
    eAdjusted = L6Convert(E, fast.EKp);
    scale = eAdjusted ./ E;
    Suse = S .* scale;
    Cuse = C .* scale;
    Iuse = InhMulp(I, fast.IKp);
    Euse = 0.6923 .* Suse + 0.3077 .* Cuse;

    recurrentE = fast.AE * [Suse; Cuse];
    recurrentI = fast.AI * Iuse;

    field = reshape(Euse, fast.fieldRows, fast.fieldCols);
    padded = padarray(field, [1 1], 'circular');
    filtered = conv2(padded, fast.kernel, 'same');
    l6Input = filtered(2:end-1, 2:end-1);
    if isempty(fast.l6PP)
        l6Axis = L6Convert(l6Input, fast.l6Pars);
    else
        l6Axis = ppval(fast.l6PP, l6Input);
    end
    l6Axis = min(max(l6Axis(:), 3), fast.l6Max) ./ 3;

    dynamicL6 = cast(repmat(l6Axis, 5, 1), fast.nnPrecision);
    dynamicE = cast(repmat(high_l4e_guard(reshape(recurrentE, n, 3)), 5, 1), fast.nnPrecision);
    dynamicI = cast(repmat(reshape(recurrentI, n, 3), 5, 1), fast.nnPrecision);
    if strcmp(fast.nnPrecision, 'single')
        output = h96_nn_predict_pages_single_mex(dynamicL6, dynamicE, dynamicI, ...
            fast.staticPages, fast.mu, fast.sd, fast.W1T, fast.b1, fast.W2T, fast.b2, ...
            fast.W3T, fast.b3, fast.W4T, fast.b4, fast.pixelWeights(:,:,1));
    else
        output = h96_nn_predict_pages_double_mex(dynamicL6, dynamicE, dynamicI, ...
            fast.staticPages, fast.mu, fast.sd, fast.W1T, fast.b1, fast.W2T, fast.b2, ...
            fast.W3T, fast.b3, fast.W4T, fast.b4, fast.pixelWeights(:,:,1));
    end
    output = double(output);
    S = fast.p .* output(:,1) + fast.oneMinusP .* S;
    C = fast.p .* output(:,2) + fast.oneMinusP .* C;
    I = fast.p .* output(:,3) + fast.oneMinusP .* I;
end

info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds / iterations;
state = struct('S', S, 'C', C, 'I', I);
end

function xg = high_l4e_guard(x)
edge = 47500;
full = 55000;
tailSlope = 0.02;
u = min(max((x - edge) ./ (full - edge), 0), 1);
g = u.^3 .* (10 - 15 .* u + 6 .* u.^2);
xg = (1 - g) .* x + g .* (edge + tailSlope .* (x - edge));
end
