function inspect_spatial_kernel_ranges(setupFile)
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
names={'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II'};
for index=1:numel(names)
    matrix=sparse(context.(names{index}));
    row=full(matrix(1,:));
    [y,x]=ind2sub(context.MapSize,1:numel(row));
    dx=mod(x-1+context.MapSize(2)/2,context.MapSize(2))-context.MapSize(2)/2;
    dy=mod(y-1+context.MapSize(1)/2,context.MapSize(1))-context.MapSize(1)/2;
    mass=sum(abs(row));
    rmsRadius=sqrt(sum(abs(row).*(dx.^2+dy.^2))/mass);
    nonzero=row~=0;
    fprintf('%s nnz/row=%d rowSum=%.12g rmsRadius=%.6f max=%.6g minNonzero=%.6g maxRadius=%.6f\n', ...
        names{index},nnz(row),sum(row),rmsRadius,max(row), ...
        min(row(nonzero)),max(sqrt(dx(nonzero).^2+dy(nonzero).^2)));
end
fprintf('L6 kernel size %s sum %.12g\n', ...
    mat2str(size(context.L6Kernel)),sum(context.L6Kernel(:)));
end
