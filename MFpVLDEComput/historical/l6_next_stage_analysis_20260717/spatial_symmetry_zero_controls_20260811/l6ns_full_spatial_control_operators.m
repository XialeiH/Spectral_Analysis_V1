function [operators,metadata] = l6ns_full_spatial_control_operators(context,control)
% Build one complete nonlinear spatial-connectivity intervention.

names={'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II'};
for index=1:numel(names)
    operators.(names{index})=sparse(context.(names{index}));
end
operators.L6=build_conv_matrix_circular(context.L6Kernel, ...
    context.MapSize(1),context.MapSize(2));
metadata=struct('Control',string(control),'L4Alpha',0,'L6Alpha',0, ...
    'PermutationSeed',NaN,'Definition',"baseline");

switch string(control)
    case "baseline"
    case "shared_source_permutation"
        rng(8112026,'twister');
        order=randperm(prod(context.MapSize));
        for index=1:numel(names)
            operators.(names{index})=operators.(names{index})(:,order);
        end
        operators.L6=operators.L6(:,order);
        metadata.PermutationSeed=8112026;
        metadata.Definition="same source-pixel permutation in every L4 block and L6";
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
    case "L4_left_half"
        [shapeOperators,shapeDetails]=l6ns_spatial_shape_operators( ...
            context,'left_half','baseline');
        for index=1:numel(names); operators.(names{index})=shapeOperators.(names{index}); end
        metadata.ShapeDetails=shapeDetails;
        metadata.Definition="one-sided half-Gaussian L4 footprints; row sums preserved";
    case "L4_triangle"
        [shapeOperators,shapeDetails]=l6ns_spatial_shape_operators( ...
            context,'triangle','baseline');
        for index=1:numel(names); operators.(names{index})=shapeOperators.(names{index}); end
        metadata.ShapeDetails=shapeDetails;
        metadata.Definition="directed triangular L4 footprints; row sums preserved";
    otherwise
        error('SpatialControl:Unknown','Unknown control %s.',control);
end

metadata.MaximumRowSumRelativeError=0;
for index=1:numel(names)
    original=context.(names{index});
    errorValue=norm(sum(operators.(names{index}),2)-sum(original,2))/ ...
        max(norm(sum(original,2)),eps);
    metadata.MaximumRowSumRelativeError=max( ...
        metadata.MaximumRowSumRelativeError,errorValue);
end
originalK=build_conv_matrix_circular(context.L6Kernel, ...
    context.MapSize(1),context.MapSize(2));
metadata.L6RowSumRelativeError=norm(sum(operators.L6,2)-sum(originalK,2))/ ...
    max(norm(sum(originalK,2)),eps);
end

function [operators,metadata]=local_flatten(operators,names,l4Alpha,l6Alpha,metadata)
n=size(operators.L6,1);
for index=1:numel(names)
    matrix=operators.(names{index});
    rowSums=full(sum(matrix,2));
    if max(rowSums)-min(rowSums)>1e-10*max(1,max(abs(rowSums)))
        error('SpatialControl:RowSums','Expected translation-invariant row sums.');
    end
    operators.(names{index})=(1-l4Alpha)*matrix+ ...
        l4Alpha*mean(rowSums)*speye(n);
end
rowSums=full(sum(operators.L6,2));
operators.L6=(1-l6Alpha)*operators.L6+l6Alpha*mean(rowSums)*speye(n);
metadata.L4Alpha=l4Alpha;
metadata.L6Alpha=l6Alpha;
metadata.Definition=sprintf('L4 alpha %.2f; L6 alpha %.2f; row sums preserved', ...
    l4Alpha,l6Alpha);
end
