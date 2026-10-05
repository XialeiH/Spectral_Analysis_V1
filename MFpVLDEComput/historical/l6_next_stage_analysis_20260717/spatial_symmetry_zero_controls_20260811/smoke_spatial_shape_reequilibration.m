function smoke_spatial_shape_reequilibration(setupFile)
% Verify controlled matrices enter the actual h96 iteration safely.

loaded=load(setupFile,'setup');
context=loaded.setup.Context;
[baseline,baselineDetails]=l6ns_spatial_shape_operators(context,'baseline','baseline');
names={'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II'};
matrixError=0;
for index=1:numel(names)
    matrixError=max(matrixError,norm(baseline.(names{index})- ...
        context.(names{index}),'fro')/max(norm(context.(names{index}),'fro'),eps));
end
assert(matrixError<1e-12,'Baseline reconstruction error %.3e.',matrixError);
assert(baselineDetails.L6RowSumRelativeError<1e-12);

for shape=["left_half","triangle"]
    [operators,details]=l6ns_spatial_shape_operators(context,shape,'baseline');
    assert(details.L4RowSumRelativeError<1e-12, ...
        '%s row-sum error %.3e.',shape,details.L4RowSumRelativeError);
    state=local_unpack(context.FixedPoint(:),prod(context.MapSize));
    placeholders={[]};
    emptyLibrary={struct('S',{{}},'C',{{}},'I',{{}})};
    [history,~,~,~,~,~,nanFlag]= ...
        LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle( ...
        context.PixLGNCtgr,context.L6Kernel,state,0.33, ...
        context.L6Parameters,5, ...
        operators.C_SS,operators.C_CS,operators.C_IS, ...
        operators.C_SC,operators.C_CC,operators.C_IC, ...
        operators.C_SI,operators.C_CI,operators.C_II, ...
        context.L4SEp,context.L4SIp,context.L4CEp,context.L4CIp, ...
        context.L4IEp,context.L4IIp,placeholders,placeholders,emptyLibrary, ...
        context.ContrastUse,context.OrientationUse,4,10,10, ...
        context.Isaturation,'xn',context.EKpUse,context.IKpUse);
    final=local_pack(history{end});
    assert(~nanFlag && all(isfinite(final)),'%s iteration failed.',shape);
    fprintf('%s: row-sum %.3e, five-step relative movement %.6e\n', ...
        shape,details.L4RowSumRelativeError, ...
        norm(final-context.FixedPoint(:))/norm(context.FixedPoint(:)));
end
fprintf('SPATIAL_REEQUILIBRATION_SMOKE_OK baseline reconstruction %.3e\n',matrixError);
end

function state=local_unpack(vector,n)
state=struct('S',vector(1:n),'C',vector(n+(1:n)),'I',vector(2*n+(1:n)));
end

function vector=local_pack(state)
vector=[state.S(:);state.C(:);state.I(:)];
end
