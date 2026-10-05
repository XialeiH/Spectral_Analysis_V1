function [state, info] = h96_nn_fixedpoint_optimized(initialState, fast, iterations, directCircular)
% Page-batched iteration with precompiled E/I nonlinearities.
if nargin<4
    directCircular = false;
end
S = initialState.S;
C = initialState.C;
I = initialState.I;
n = fast.n;

timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S+0.3077.*C;
    eAdjusted = h96_e_mulp_fast(E, fast.ekpFast);
    scale = eAdjusted./E;
    Suse = S.*scale;
    Cuse = C.*scale;
    Iuse = h96_inh_mulp_fast(I, fast.ikpFast);
    Euse = 0.6923.*Suse+0.3077.*Cuse;

    recurrentE = fast.AE*[Suse; Cuse];
    recurrentI = fast.AI*Iuse;

    field = reshape(Euse, fast.fieldRows, fast.fieldCols);
    if directCircular
        l6Input = h96_circular_conv3x3(field, fast.kernel);
    else
        padded = padarray(field, [1 1], 'circular');
        filtered = conv2(padded, fast.kernel, 'same');
        l6Input = filtered(2:end-1, 2:end-1);
    end
    if isempty(fast.l6PP)
        l6Axis = L6Convert(l6Input, fast.l6Pars);
    else
        l6Axis = ppval(fast.l6PP, l6Input);
    end
    l6Axis = min(max(l6Axis(:), 3), fast.l6Max)./3;

    output = h96_nn_all_population_responses(recurrentE, recurrentI, l6Axis, fast);
    S = fast.p.*output(:,1)+fast.oneMinusP.*S;
    C = fast.p.*output(:,2)+fast.oneMinusP.*C;
    I = fast.p.*output(:,3)+fast.oneMinusP.*I;
end

info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds/iterations;
state = struct('S', S, 'C', C, 'I', I);
end
