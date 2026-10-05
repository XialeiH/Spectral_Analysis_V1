function [state, info] = h96_nn_fixedpoint_compact(initialState, fast, iterations)
% NN-only iteration with compact nonzero LGN pairs and page-batched S/C/I.
S = initialState.S;
C = initialState.C;
I = initialState.I;

timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S+0.3077.*C;
    eAdjusted = L6Convert(E, fast.EKp);
    scale = eAdjusted./E;
    Suse = S.*scale;
    Cuse = C.*scale;
    Iuse = InhMulp(I, fast.IKp);
    Euse = 0.6923.*Suse+0.3077.*Cuse;

    recurrentE = fast.AE*[Suse; Cuse];
    recurrentI = fast.AI*Iuse;

    field = reshape(Euse, fast.fieldRows, fast.fieldCols);
    padded = padarray(field, [1 1], 'circular');
    filtered = conv2(padded, fast.kernel, 'same');
    l6Input = filtered(2:end-1, 2:end-1);
    if isempty(fast.l6PP)
        l6Axis = L6Convert(l6Input, fast.l6Pars);
    else
        l6Axis = ppval(fast.l6PP, l6Input);
    end
    l6Axis = min(max(l6Axis(:), 3), fast.l6Max)./3;

    output = h96_nn_all_population_responses_compact(recurrentE, recurrentI, l6Axis, fast);
    S = fast.p.*output(:,1)+fast.oneMinusP.*S;
    C = fast.p.*output(:,2)+fast.oneMinusP.*C;
    I = fast.p.*output(:,3)+fast.oneMinusP.*I;
end

info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds/iterations;
state = struct('S', S, 'C', C, 'I', I);
end
