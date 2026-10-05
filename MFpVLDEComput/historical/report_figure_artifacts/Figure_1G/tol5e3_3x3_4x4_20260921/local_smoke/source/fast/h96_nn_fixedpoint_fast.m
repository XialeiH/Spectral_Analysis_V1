function [state, info] = h96_nn_fixedpoint_fast(initialState, fast, iterations)
% NN-only fixed-count iteration with no library state or trace allocation.
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

    Sout = population_response(recurrentE(1:n), recurrentI(1:n), l6Axis, fast, 'S');
    Cout = population_response(recurrentE(n+1:2*n), recurrentI(n+1:2*n), l6Axis, fast, 'C');
    Iout = population_response(recurrentE(2*n+1:3*n), recurrentI(2*n+1:3*n), l6Axis, fast, 'I');

    S = fast.p .* Sout + fast.oneMinusP .* S;
    C = fast.p .* Cout + fast.oneMinusP .* C;
    I = fast.p .* Iout + fast.oneMinusP .* I;
end

info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds / iterations;
state = struct('S', S, 'C', C, 'I', I);
end

function output = population_response(l4e, l4i, l6Axis, fast, population)
model = fast.models.(population);
x = fast.staticNormalized.(population);
x(:,3) = (repmat(l6Axis, 5, 1) - model.mu(3)) / model.sd(3);
x(:,4) = (repmat(high_l4e_guard(l4e), 5, 1) - model.mu(4)) / model.sd(4);
x(:,5) = (repmat(l4i, 5, 1) - model.mu(5)) / model.sd(5);
y = reshape(h96_nn_predict_cached(x, model), fast.n, 5);
output = sum(fast.PixLGNCtgr .* y, 2);
end

function xg = high_l4e_guard(x)
edge = 47500;
full = 55000;
tailSlope = 0.02;
u = min(max((x - edge) ./ (full - edge), 0), 1);
g = u.^3 .* (10 - 15 .* u + 6 .* u.^2);
xg = (1 - g) .* x + g .* (edge + tailSlope .* (x - edge));
end
