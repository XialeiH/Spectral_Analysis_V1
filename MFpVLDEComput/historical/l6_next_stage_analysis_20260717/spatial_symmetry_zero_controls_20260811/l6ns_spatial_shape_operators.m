function [operators,details] = l6ns_spatial_shape_operators(context,l4Shape,l6Shape)
% Construct row-sum-preserving spatial-footprint controls.

l4Names = {'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II'};
l4Controlled = cell(size(l4Names));
l4Details = cell(size(l4Names));
for index = 1:numel(l4Names)
    [l4Controlled{index},l4Details{index}] = local_shape_operator( ...
        context.(l4Names{index}),l4Shape,context.MapSize);
end

kBase = build_conv_matrix_circular(context.L6Kernel, ...
    context.MapSize(1),context.MapSize(2));
[kControl,l6Details] = local_shape_operator(kBase,l6Shape,context.MapSize);

operators = cell2struct(l4Controlled,l4Names,2);
operators.L6 = kControl;

details = struct();
details.L4Shape = string(l4Shape);
details.L6Shape = string(l6Shape);
details.L4RowSumRelativeError = max(cellfun( ...
    @(item) item.RowSumRelativeError,l4Details));
details.L4TranslationError = mean(cellfun( ...
    @(item) item.TranslationError,l4Details));
details.L4Rotation90Error = mean(cellfun( ...
    @(item) item.Rotation90Error,l4Details));
details.L4Rotation180Error = mean(cellfun( ...
    @(item) item.Rotation180Error,l4Details));
details.L4ReflectionError = mean(cellfun( ...
    @(item) item.ReflectionError,l4Details));
details.L4RetainedMassFraction = mean(cellfun( ...
    @(item) item.RetainedMassFraction,l4Details));
details.L4TriangleRadius = mean(cellfun( ...
    @(item) item.TriangleRadius,l4Details));
details.L6RowSumRelativeError = l6Details.RowSumRelativeError;
details.L6TranslationError = l6Details.TranslationError;
details.L6Rotation90Error = l6Details.Rotation90Error;
details.L6Rotation180Error = l6Details.Rotation180Error;
details.L6ReflectionError = l6Details.ReflectionError;
details.L6RetainedMassFraction = l6Details.RetainedMassFraction;
details.L6TriangleRadius = l6Details.TriangleRadius;
details.L4Kernel = local_centered_kernel(l4Controlled{1},context.MapSize);
details.L6Kernel = local_centered_kernel(kControl,context.MapSize);
end

function [controlled,details] = local_shape_operator(matrix,shape,mapSize)
matrix = sparse(matrix);
n = size(matrix,1);
[rows,columns,values] = find(matrix);
[rowY,rowX] = ind2sub(mapSize,rows);
[columnY,columnX] = ind2sub(mapSize,columns);
dx = local_periodic_offset(columnX-rowX,mapSize(2));
dy = local_periodic_offset(columnY-rowY,mapSize(1));

switch string(shape)
    case "baseline"
        keep = true(size(values));
        triangleRadius = NaN;
    case "left_half"
        keep = dx<0 | (dx==0 & dy==0);
        triangleRadius = NaN;
    case "triangle"
        leftKeep = dx<0 | (dx==0 & dy==0);
        targetFraction = sum(abs(values(leftKeep)))/sum(abs(values));
        radiusCandidates = 0:0.25:(2*max(mapSize));
        retained = zeros(size(radiusCandidates));
        for radiusIndex = 1:numel(radiusCandidates)
            candidateKeep = local_triangle_mask(dx,dy,radiusCandidates(radiusIndex));
            retained(radiusIndex) = sum(abs(values(candidateKeep)))/sum(abs(values));
        end
        [~,bestIndex] = min(abs(retained-targetFraction));
        triangleRadius = radiusCandidates(bestIndex);
        keep = local_triangle_mask(dx,dy,triangleRadius);
    otherwise
        error('SpatialShape:Unknown','Unknown spatial shape %s.',shape);
end

masked = sparse(rows(keep),columns(keep),values(keep),n,n);
targetRowSums = full(sum(matrix,2));
maskedRowSums = full(sum(masked,2));
if any(abs(maskedRowSums)<eps & abs(targetRowSums)>eps)
    error('SpatialShape:EmptyRow','Shape %s removed an entire nonzero row.',shape);
end
rowScale = targetRowSums./maskedRowSums;
rowScale(~isfinite(rowScale)) = 1;
controlled = spdiags(rowScale,0,n,n)*masked;

details = local_operator_details(matrix,controlled,mapSize);
details.RetainedMassFraction = sum(abs(values(keep)))/sum(abs(values));
details.TriangleRadius = triangleRadius;
end

function keep = local_triangle_mask(dx,dy,radius)
% Centered equilateral triangle pointing right.
keep = dx>=-radius/2 & dx<=radius & abs(dy)<=(radius-dx)/sqrt(3);
end

function offset = local_periodic_offset(offset,period)
offset = mod(offset+floor(period/2),period)-floor(period/2);
end

function details = local_operator_details(original,controlled,mapSize)
n = size(controlled,1);
grid = reshape(1:n,mapSize);
translationOrder = circshift(grid,[1 0]); translationOrder = translationOrder(:);
rotation90Order = rot90(grid,1); rotation90Order = rotation90Order(:);
rotation180Order = rot90(grid,2); rotation180Order = rotation180Order(:);
reflectionOrder = fliplr(grid); reflectionOrder = reflectionOrder(:);
denominator = max(norm(controlled,'fro'),eps);
details.RowSumRelativeError = norm(sum(controlled,2)-sum(original,2))/ ...
    max(norm(sum(original,2)),eps);
details.TranslationError = norm(controlled(translationOrder,translationOrder)- ...
    controlled,'fro')/denominator;
details.Rotation90Error = norm(controlled(rotation90Order,rotation90Order)- ...
    controlled,'fro')/denominator;
details.Rotation180Error = norm(controlled(rotation180Order,rotation180Order)- ...
    controlled,'fro')/denominator;
details.ReflectionError = norm(controlled(reflectionOrder,reflectionOrder)- ...
    controlled,'fro')/denominator;
end

function kernel = local_centered_kernel(matrix,mapSize)
kernel = reshape(full(matrix(1,:)),mapSize);
kernel = circshift(kernel,[floor(mapSize(1)/2),floor(mapSize(2)/2)]);
kernel = kernel/max(sum(abs(kernel(:))),eps);
end
