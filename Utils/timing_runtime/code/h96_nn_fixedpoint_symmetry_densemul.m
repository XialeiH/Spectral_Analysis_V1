function info = h96_nn_fixedpoint_symmetry_densemul(initialState,fast,iterations)
% h96 iteration with dense-block recurrent products and no sparse multiply.
S = initialState.S;
C = initialState.C;
I = initialState.I;
timer = tic;
for epoch = 1:iterations
    E = 0.6923.*S+0.3077.*C;
    adjustedE = L6Convert(E,fast.EKp);
    scale = adjustedE./E;
    Suse = S.*scale;
    Cuse = C.*scale;
    Iuse = InhMulp(I,fast.IKp);
    Euse = 0.6923.*Suse+0.3077.*Cuse;
    recurrentE = dense_matvec(fast.AE,[Suse;Cuse]);
    recurrentI = dense_matvec(fast.AI,Iuse);
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
    output = h96_nn_all_population_responses_symmetry_nosparse( ...
        recurrentE,recurrentI,l6Axis,fast);
    S = fast.p.*output(:,1)+fast.oneMinusP.*S;
    C = fast.p.*output(:,2)+fast.oneMinusP.*C;
    I = fast.p.*output(:,3)+fast.oneMinusP.*I;
end
info.Iterations = iterations;
info.Seconds = toc(timer);
info.SecondsPerIteration = info.Seconds/iterations;
info.SymmetryClassCount = fast.symClassCount;
info.Checksum = sum(S)+sum(C)+sum(I);
end

function y = dense_matvec(A,x)
nRows = size(A,1);
nCols = size(A,2);
targetBytes = 256*1024^2;
rowsPerBlock = max(1,floor(targetBytes/(8*nCols)));
y = zeros(nRows,1);
for first = 1:rowsPerBlock:nRows
    rows = first:min(first+rowsPerBlock-1,nRows);
    denseBlock = full(A(rows,:));
    y(rows) = denseBlock*x;
end
end
