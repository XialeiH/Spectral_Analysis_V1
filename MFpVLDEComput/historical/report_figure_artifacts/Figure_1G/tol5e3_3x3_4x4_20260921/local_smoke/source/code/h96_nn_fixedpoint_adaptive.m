function [state,info] = h96_nn_fixedpoint_adaptive(initialState,fast,iterations)
% Discover exact repeated tuples once, then reuse those spatial classes.
S = initialState.S;
C = initialState.C;
I = initialState.I;
cache = [];
timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S+0.3077.*C;
    adjustedE = L6Convert(E,fast.EKp);
    scale = adjustedE./E;
    Suse = S.*scale;
    Cuse = C.*scale;
    Iuse = InhMulp(I,fast.IKp);
    Euse = 0.6923.*Suse+0.3077.*Cuse;
    recurrentE = fast.AE*[Suse;Cuse];
    recurrentI = fast.AI*Iuse;
    field = reshape(Euse,fast.fieldRows,fast.fieldCols);
    padded = padarray(field,[1 1],'circular');
    filtered = conv2(padded,fast.kernel,'same');
    l6Input = filtered(2:end-1,2:end-1);
    if isempty(fast.l6PP)
        l6Axis = L6Convert(l6Input,fast.l6Pars);
    else
        l6Axis = ppval(fast.l6PP,l6Input);
    end
    l6Axis = min(max(l6Axis(:),3),fast.l6Max)./3;
    if epoch==1
        [output,uniqueCounts,cache] = h96_nn_all_population_responses_unique( ...
            recurrentE,recurrentI,l6Axis,fast);
    else
        output = h96_nn_all_population_responses_cached( ...
            recurrentE,recurrentI,l6Axis,fast,cache);
    end
    S = fast.p.*output(:,1)+fast.oneMinusP.*S;
    C = fast.p.*output(:,2)+fast.oneMinusP.*C;
    I = fast.p.*output(:,3)+fast.oneMinusP.*I;
end
info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds/iterations;
info.UniqueCounts = uniqueCounts;
state = struct('S',S,'C',C,'I',I);
end
