function curves = figure6_0_tuning_curves(canonicalStates,canonicalAngles,cWeight)
% Reconstruct the two project E tuning curves from seven canonical states.

fullAngles = 0:3.75:180;
pixelRows = [5 1];
pixelCols = [10 10];
mapSide = 40;
n = mapSide^2;
curves = zeros(numel(fullAngles),2);
for angleIndex = 1:numel(fullAngles)
    angle = fullAngles(angleIndex);
    for pixelIndex = 1:2
        [canonicalAngle,row,col] = local_map_angle_pixel( ...
            angle,pixelRows(pixelIndex),pixelCols(pixelIndex));
        canonicalIndex = find(abs(canonicalAngles-canonicalAngle)<1e-10,1);
        if isempty(canonicalIndex)
            error('Figure60:TuningAngle','Missing canonical angle %.6g.',canonicalAngle);
        end
        state = canonicalStates(:,canonicalIndex);
        sMap = reshape(state(1:n),mapSide,mapSide);
        cMap = reshape(state(n+(1:n)),mapSide,mapSide);
        curves(angleIndex,pixelIndex) = ...
            (1-cWeight)*sMap(row,col)+cWeight*cMap(row,col);
    end
end
end

function [canonicalAngle,row,col] = local_map_angle_pixel(angle,row,col)
if row==5 && col==10
    if angle<=22.5
        canonicalAngle=angle; row=5; col=10;
    elseif angle<45
        canonicalAngle=45-angle; row=1; col=6;
    elseif angle==45
        canonicalAngle=0; row=10; col=6;
    elseif angle<=67.5
        canonicalAngle=angle-45; row=10; col=6;
    elseif angle<=90
        canonicalAngle=90-angle; row=5; col=1;
    elseif angle<112.5
        canonicalAngle=angle-90; row=5; col=1;
    elseif angle<135
        canonicalAngle=135-angle; row=10; col=6;
    elseif angle<157.5
        canonicalAngle=angle-135; row=1; col=6;
    else
        canonicalAngle=180-angle; row=5; col=10;
    end
    return
end
if row==1 && col==10
    if angle<=22.5
        canonicalAngle=angle; row=1; col=10;
    elseif angle<=45
        canonicalAngle=45-angle; row=1; col=10;
    elseif angle<=67.5
        canonicalAngle=angle-45; row=10; col=10;
    elseif angle<=90
        canonicalAngle=90-angle; row=1; col=1;
    elseif angle<=112.5
        canonicalAngle=angle-90; row=10; col=1;
    elseif angle<=135
        canonicalAngle=135-angle; row=10; col=1;
    elseif angle<=157.5
        canonicalAngle=angle-135; row=1; col=1;
    else
        canonicalAngle=180-angle; row=10; col=10;
    end
    return
end
error('Figure60:TuningPixel','Unsupported tuning pixel.');
end
