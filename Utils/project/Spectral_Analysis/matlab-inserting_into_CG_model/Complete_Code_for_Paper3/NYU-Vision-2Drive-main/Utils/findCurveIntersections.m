function intersections = findCurveIntersections(x, y1, y2)
    % Ensure x, y1, and y2 are column vectors
    x = x(:);
    y1 = y1(:);
    y2 = y2(:);

    % Calculate the differences between the two y-data vectors
    diffY = y1 - y2;

    % Find indices where the difference changes sign
    signChangeIndices = find(diff(diffY > 0));

    % Preallocate the intersections array
    intersections = zeros(length(signChangeIndices), 2);

    % Loop through each sign change index to compute the intersection points
    for i = 1:length(signChangeIndices)
        % Linear interpolation to find a more accurate x value for the intersection
        idx = signChangeIndices(i);
        x1 = x(idx);
        x2 = x(idx + 1);
        y1_1 = y1(idx);
        y1_2 = y1(idx + 1);
        y2_1 = y2(idx);
        y2_2 = y2(idx + 1);

        % Linear interpolation formula
        xIntersect = x1 + (x2 - x1) * abs(y1_1 - y2_1) / abs((y1_1 - y2_1) - (y1_2 - y2_2));
        
        % Compute corresponding y value for the intersection point
        yIntersect = interp1([x1, x2], [y1_1, y1_2], xIntersect, 'linear');

        % Store the computed intersection
        intersections(i, :) = [xIntersect, yIntersect];
    end

    % Remove any rows with NaNs, if any
    intersections = intersections(~any(isnan(intersections), 2), :);
end