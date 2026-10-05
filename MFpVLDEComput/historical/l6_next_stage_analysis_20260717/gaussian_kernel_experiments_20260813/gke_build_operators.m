function [operators,metadata,audit]=gke_build_operators(context,controlName)
% Build one controlled connectivity model while isolating the named change.

names={'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II'};
for index=1:numel(names)
    baseline.(names{index})=sparse(context.(names{index}));
end
baseline.L6=build_conv_matrix_circular(context.L6Kernel, ...
    context.MapSize(1),context.MapSize(2));
operators=baseline;
controlName=string(controlName);
metadata=struct('Control',controlName,'Seed',NaN,'Parameter',NaN, ...
    'Definition',"baseline",'ExpectedRowSumPreservation',true, ...
    'MaskRetainedMassFraction',NaN,'StrengthFactorCV',0, ...
    'CWeight',context.CWeight);

switch controlName
    case "baseline"
    case "swap_EI_ranges"
        operators=local_swap_ei_ranges(operators,context.CWeight);
        metadata.Definition=["swap normalized E-source and I-source spatial " ...
            "profiles within each postsynaptic population; preserve each block row sum"];
    case "shape_inscribed_triangle"
        [operators,metadata]=local_apply_shape(operators,names,"triangle",metadata,context.MapSize);
    case "shape_inscribed_square"
        [operators,metadata]=local_apply_shape(operators,names,"square",metadata,context.MapSize);
    case "shape_inscribed_hexagon"
        [operators,metadata]=local_apply_shape(operators,names,"hexagon",metadata,context.MapSize);
    case "shape_left_half"
        [operators,metadata]=local_apply_shape(operators,names,"left_half",metadata,context.MapSize);
    case "joint_alpha_0p10"
        [operators,metadata]=local_flatten(operators,names,0.10,0.10,metadata);
    case "joint_alpha_0p25"
        [operators,metadata]=local_flatten(operators,names,0.25,0.25,metadata);
    case "joint_alpha_0p50"
        [operators,metadata]=local_flatten(operators,names,0.50,0.50,metadata);
    case "remove_L6_smoothing"
        [operators,metadata]=local_flatten(operators,names,0,1,metadata);
    case "remove_L4_smoothing"
        [operators,metadata]=local_flatten(operators,names,1,0,metadata);
    case "remove_L4_and_L6_smoothing"
        [operators,metadata]=local_flatten(operators,names,1,1,metadata);
    case "translation_jitter_1px"
        [operators,metadata]=local_target_jitter(operators,names,1,813101,metadata,context.MapSize);
    case "translation_jitter_2px_mix0p50"
        original=operators;
        [shifted,metadata]=local_target_jitter(operators,names,2,813103,metadata,context.MapSize);
        for index=1:numel(names)
            operators.(names{index})=0.5*original.(names{index})+ ...
                0.5*shifted.(names{index});
        end
        metadata.Parameter=[2 0.5];
        metadata.Definition=["50% mixture with targetwise random kernel-center " ...
            "jitter up to 2 pixels; shared displacement across L4 blocks; " ...
            "row sums preserved"];
    case "translation_sinusoidal_warp_2px"
        [operators,metadata]=local_sinusoidal_warp(operators,names,2,metadata,context.MapSize);
    case "translation_source_permutation_0p25"
        rng(813116,'twister');
        order=randperm(prod(context.MapSize));
        for index=1:numel(names)
            matrix=operators.(names{index});
            operators.(names{index})=0.75*matrix+0.25*matrix(:,order);
        end
        metadata.Seed=813116;
        metadata.Parameter=0.25;
        metadata.Definition=["25% mixture with the same random source-pixel " ...
            "permutation in all L4 blocks; " ...
            "row sums preserved, locality and translation invariance destroyed"];
    case "kernel_noise_cv0p10"
        [operators,metadata]=local_kernel_noise(operators,names,0.10,813171,metadata,context.MapSize);
    case "kernel_noise_cv0p30"
        [operators,metadata]=local_kernel_noise(operators,names,0.30,813173,metadata,context.MapSize);
    case "strength_noise_cv0p05"
        [operators,metadata]=local_strength_noise(operators,names,0.05,813191,metadata);
    case "strength_noise_cv0p15"
        [operators,metadata]=local_strength_noise(operators,names,0.15,813193,metadata);
    otherwise
        error('GKE:UnknownControl','Unknown control %s.',controlName);
end

audit=gke_operator_audit(baseline,operators,context.MapSize);
metadata.MaximumL4RowSumRelativeError=max(audit.rowSumRelativeError(1:9));
metadata.L6RowSumRelativeError=audit.rowSumRelativeError(10);
metadata.MeanL4RowSumRatio=mean(audit.meanRowSumRatio(1:9));
metadata.MeanL4RowSumCV=mean(audit.controlledRowSumCV(1:9));
metadata.MeanL4TranslationError=mean(audit.translationError(1:9));
metadata.MeanL4Rotation90Error=mean(audit.rotation90Error(1:9));
end

function operators=local_swap_ei_ranges(operators,cWeight)
triples={{'C_SS','C_SC','C_SI'},{'C_CS','C_CC','C_CI'}, ...
    {'C_IS','C_IC','C_II'}};
for target=1:numel(triples)
    names=triples{target};
    eS=operators.(names{1}); eC=operators.(names{2}); inhibitory=operators.(names{3});
    eProfile=(1-cWeight)*local_row_normalize(eS)+cWeight*local_row_normalize(eC);
    iProfile=local_row_normalize(inhibitory);
    operators.(names{1})=local_apply_profile(iProfile,eS);
    operators.(names{2})=local_apply_profile(iProfile,eC);
    operators.(names{3})=local_apply_profile(eProfile,inhibitory);
end
end

function normalized=local_row_normalize(matrix)
rowSums=full(sum(matrix,2));
normalized=spdiags(1./rowSums,0,size(matrix,1),size(matrix,1))*matrix;
end

function result=local_apply_profile(profile,target)
result=spdiags(full(sum(target,2)),0,size(target,1),size(target,1))*profile;
end

function [operators,metadata]=local_apply_shape(operators,names,shape,metadata,mapSize)
retained=zeros(numel(names),1);
for index=1:numel(names)
    [operators.(names{index}),retained(index)]= ...
        local_mask_shape(operators.(names{index}),shape,mapSize);
end
metadata.MaskRetainedMassFraction=mean(retained);
metadata.Definition=sprintf(['keep original L4 weights inside the inscribed %s; ' ...
    'remove weights outside and rescale every row to its original sum'],shape);
end

function [controlled,retainedFraction]=local_mask_shape(matrix,shape,mapSize)
[rows,columns,values]=find(matrix);
[rowY,rowX]=ind2sub(mapSize,rows);
[columnY,columnX]=ind2sub(mapSize,columns);
dx=local_periodic_offset(columnX-rowX,mapSize(2));
dy=local_periodic_offset(columnY-rowY,mapSize(1));
radius=max(sqrt(double(dx).^2+double(dy).^2));
switch string(shape)
    case "left_half"
        keep=dx<=0;
    otherwise
        switch string(shape)
            case "triangle"; sideCount=3;
            case "square"; sideCount=4;
            case "hexagon"; sideCount=6;
        end
        angles=(0:sideCount-1)*(2*pi/sideCount);
        vertexX=radius*cos(angles);
        vertexY=radius*sin(angles);
        [inside,onBoundary]=inpolygon(double(dx),double(dy),vertexX,vertexY);
        keep=inside|onBoundary;
end
retainedFraction=sum(abs(values(keep)))/sum(abs(values));
controlled=local_rescaled_sparse(rows(keep),columns(keep),values(keep),matrix);
end

function [operators,metadata]=local_flatten(operators,names,l4Alpha,l6Alpha,metadata)
n=size(operators.L6,1);
for index=1:numel(names)
    matrix=operators.(names{index});
    rowSums=full(sum(matrix,2));
    operators.(names{index})=(1-l4Alpha)*matrix+ ...
        l4Alpha*mean(rowSums)*speye(n);
end
rowSums=full(sum(operators.L6,2));
operators.L6=(1-l6Alpha)*operators.L6+l6Alpha*mean(rowSums)*speye(n);
metadata.Parameter=max(l4Alpha,l6Alpha);
metadata.Definition=sprintf('L4 alpha %.2f, L6 alpha %.2f; row sums preserved', ...
    l4Alpha,l6Alpha);
end

function [operators,metadata]=local_target_jitter(operators,names,maxShift,seed,metadata,mapSize)
rng(seed,'twister');
n=prod(mapSize);
shiftY=randi([-maxShift maxShift],n,1);
shiftX=randi([-maxShift maxShift],n,1);
for index=1:numel(names)
    operators.(names{index})=local_shift_rows(operators.(names{index}),shiftY,shiftX,mapSize);
end
metadata.Seed=seed;
metadata.Parameter=maxShift;
metadata.Definition=sprintf(['same random targetwise center displacement across all L4 blocks; ' ...
    'integer shifts in [-%d,%d] pixels; row weights and sums preserved'],maxShift,maxShift);
end

function [operators,metadata]=local_sinusoidal_warp(operators,names,amplitude,metadata,mapSize)
n=prod(mapSize);
[targetY,targetX]=ind2sub(mapSize,(1:n)');
shiftX=round(amplitude*sin(2*pi*(targetY-1)/mapSize(1)));
shiftY=round(amplitude*sin(2*pi*(targetX-1)/mapSize(2)));
for index=1:numel(names)
    operators.(names{index})=local_shift_rows(operators.(names{index}),shiftY,shiftX,mapSize);
end
metadata.Parameter=amplitude;
metadata.Definition=sprintf(['smooth position-dependent kernel-center warp of amplitude %d pixels; ' ...
    'row weights and sums preserved'],amplitude);
end

function controlled=local_shift_rows(matrix,shiftY,shiftX,mapSize)
[rows,columns,values]=find(matrix);
[columnY,columnX]=ind2sub(mapSize,columns);
newY=mod(columnY-1+shiftY(rows),mapSize(1))+1;
newX=mod(columnX-1+shiftX(rows),mapSize(2))+1;
newColumns=sub2ind(mapSize,newY,newX);
controlled=sparse(rows,newColumns,values,size(matrix,1),size(matrix,2));
end

function [operators,metadata]=local_kernel_noise(operators,names,noiseCV,seed,metadata,mapSize)
rng(seed,'twister');
factorMap=max(0.05,1+noiseCV*randn(mapSize));
factorMap=factorMap/mean(factorMap(:));
for index=1:numel(names)
    matrix=operators.(names{index});
    [rows,columns,values]=find(matrix);
    [rowY,rowX]=ind2sub(mapSize,rows);
    [columnY,columnX]=ind2sub(mapSize,columns);
    dx=mod(columnX-rowX,mapSize(2));
    dy=mod(columnY-rowY,mapSize(1));
    factors=factorMap(sub2ind(mapSize,dy+1,dx+1));
    operators.(names{index})=local_rescaled_sparse( ...
        rows,columns,values.*factors,matrix);
end
metadata.Seed=seed;
metadata.Parameter=noiseCV;
metadata.Definition=sprintf(['one shared multiplicative offset-kernel noise template, CV %.2f; ' ...
    'translation invariance and every row sum preserved'],noiseCV);
end

function [operators,metadata]=local_strength_noise(operators,names,noiseCV,seed,metadata)
rng(seed,'twister');
n=size(operators.L6,1);
factors=max(0.10,1+noiseCV*randn(n,1));
factors=factors/mean(factors);
for index=1:numel(names)
    operators.(names{index})=spdiags(factors,0,n,n)*operators.(names{index});
end
metadata.Seed=seed;
metadata.Parameter=noiseCV;
metadata.ExpectedRowSumPreservation=false;
metadata.StrengthFactorCV=std(factors)/mean(factors);
metadata.StrengthFactors=factors;
metadata.Definition=sprintf(['shared target-pixel total-strength factors, mean 1 and CV %.2f; ' ...
    'each local kernel shape preserved, row sums intentionally heterogeneous'],noiseCV);
end

function controlled=local_rescaled_sparse(rows,columns,values,original)
n=size(original,1);
masked=sparse(rows,columns,values,n,n);
targetSums=full(sum(original,2));
maskedSums=full(sum(masked,2));
if any(abs(maskedSums)<eps & abs(targetSums)>eps)
    error('GKE:EmptyRow','Manipulation removed every connection from a nonzero row.');
end
scale=targetSums./maskedSums;
scale(~isfinite(scale))=1;
controlled=spdiags(scale,0,n,n)*masked;
end

function offset=local_periodic_offset(offset,period)
offset=mod(offset+floor(period/2),period)-floor(period/2);
end
