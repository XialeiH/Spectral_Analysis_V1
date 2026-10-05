function intersection_points = findintersection(xdata_degrees,curve_ydata,constant_y)
% Your data
%xdata_degrees = -90:7.5:90; % xdata in degrees
%curve_ydata = TuningCurvePlot; % ydata representing the tuning curve
%constant_y = 5; % Example constant y-value for the straight line

% Calculate the difference between the curve and the constant line
difference = curve_ydata - constant_y;

% Find where the sign changes in the difference array
sign_changes = diff(sign(difference));

% Find the indices where the sign changes occur
intersection_indices = find(sign_changes ~= 0);

% Initialize an array to store intersection points
intersection_points = [];

% Calculate the intersection points using linear interpolation
for idx = 1:length(intersection_indices)
    i = intersection_indices(idx);
    % Linear interpolation for x
    x1 = xdata_degrees(i);
    x2 = xdata_degrees(i+1);
    y1 = difference(i);
    y2 = difference(i+1);
    x_intersect = x1 - y1 * (x2 - x1) / (y2 - y1);
    y_intersect = constant_y; % y is constant at the intersection
    intersection_points = [intersection_points; x_intersect, y_intersect];
end

% % Display the intersection points
% disp('Intersection Points (x, y):');
% disp(intersection_points);
end