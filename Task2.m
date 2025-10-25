clc; close all;
I=imread("Object.png");
%figure, imshow(I);
IRed=double(I(:,:,1));
IGreen = double(I(:,:,2));
IBlue = double(I(:,:,3));
%take average of all three colors
IGrey = (IRed + IGreen + IBlue) / 3;
greyImage=uint8(IGrey);

%edge detection
IsobelTemp=double(greyImage);
[m,n]=size(greyImage);

Isobelx = zeros(m, n);

for i = 2:m-1
    for j = 2:n-1
        Isobelx(i, j) = (-1 * IsobelTemp(i-1, j-1)) + (0 * IsobelTemp(i-1, j)) + (1 * IsobelTemp(i-1, j+1)) ...
                      + (-2 * IsobelTemp(i, j-1))   + (0 * IsobelTemp(i, j))   + (2 * IsobelTemp(i, j+1)) ...
                      + (-1 * IsobelTemp(i+1, j-1)) + (0 * IsobelTemp(i+1, j)) + (1 * IsobelTemp(i+1, j+1));
    end
end

Isobelx = sqrt(Isobelx.^2);   % getting absolute value
Isobelx = uint8(Isobelx);


Isobely = zeros(m, n);   % y because differentiation in y direction. The edge is horizontal

for i = 2:m-1
    for j = 2:n-1
        Isobely(i, j) = ( 1 * IsobelTemp(i-1, j-1) + 2 * IsobelTemp(i-1, j) + IsobelTemp(i-1, j+1) ...
                         + 0 * IsobelTemp(i, j-1)   + 0 * IsobelTemp(i, j)   + 0 * IsobelTemp(i, j+1) ...
                         - IsobelTemp(i+1, j-1) - 2 * IsobelTemp(i+1, j) - IsobelTemp(i+1, j+1) );
    end
end

Isobely = sqrt(Isobely.^2);   % getting absolute value
Isobely = uint8(Isobely);


%combine
Isobelx = double(Isobelx);
Isobely = double(Isobely);
Isobelxy = sqrt(Isobelx.^2 + Isobely.^2);

Isobelxy = uint8(Isobelxy);
figure, imshow(Isobelxy);
title('Sobel Combined Edge Detection');

figure, imhist(Isobelxy);
title("Histogram for Thresholding")
Isobelxy=double(Isobelxy);
[m,n]=size(Isobelxy);
Ithreshold=zeros(m,n);
for i=1:m
    for j=1:n
        if Isobelxy(i,j)>40
            Ithreshold(i,j) = 255; 
        else
            Ithreshold(i,j) = 0; 
        end
    end
end

Ithreshold=uint8(Ithreshold);
% Display the thresholded image
figure, imshow(Ithreshold);
title('Thresholded Image');

I=Ithreshold;
Imedian = zeros(m, n);

for i = 2:m-1
    for j = 2:n-1
        Imedian(i, j) = median([I(i-1, j-1), I(i-1, j), I(i-1, j+1), ...
                                I(i, j-1),   I(i, j),   I(i, j+1), ...
                                I(i+1, j-1), I(i+1, j), I(i+1, j+1)]);
    end
end

Imedian = uint8(Imedian);
figure, imshow(Imedian)
title('Median Filter Applied')

% Convert thresholded image to binary (0 or 1)
binaryMask = Imedian > 0;

% Display the binary mask before drawing
figure, imshow(binaryMask);
title('Centroid and Bounding box');

hold on;
draw_box_and_centroids(binaryMask);
hold off;


figure, imshow(binaryMask);
title('Orientation and Axes');

hold on;
calculate_orientation(binaryMask);
hold off;

% Calculate the perimeter of the binary mask
calc_perimeter(binaryMask);

function calc_perimeter(binaryMask)
    [m, n] = size(binaryMask);
    
    % ---- STEP 1: Find ALL boundary pixels manually ----
    boundary = zeros(m, n);
    
    for i = 2:m-1
        for j = 2:n-1
            if binaryMask(i,j) == 1
                % Check 8-neighbors: if any neighbor is 0, this is a boundary pixel
                neighbors = [binaryMask(i-1,j-1), binaryMask(i-1,j), binaryMask(i-1,j+1), ...
                            binaryMask(i,j-1),                        binaryMask(i,j+1), ...
                            binaryMask(i+1,j-1), binaryMask(i+1,j), binaryMask(i+1,j+1)];
                
                if any(neighbors == 0)
                    boundary(i,j) = 1;
                end
            end
        end
    end
    
    % ---- STEP 2: Extract boundary coordinates ----
    [rows, cols] = find(boundary);
    
    if isempty(rows)
        disp('No boundary found');
        return;
    end
    
    % ---- STEP 3: Trace boundary in order (chain code approach) ----
    % Start from topmost-leftmost point
    [~, start_idx] = min(rows + cols/n); % prioritize top, then left
    current = [rows(start_idx), cols(start_idx)];
    
    ordered_boundary = current;
    visited = false(length(rows), 1);
    visited(start_idx) = true;
    
    % 8-connectivity search order (clockwise from top)
    directions = [-1,0; -1,1; 0,1; 1,1; 1,0; 1,-1; 0,-1; -1,-1];
    
    % Trace the boundary
    max_iterations = length(rows) * 2;
    iteration = 0;
    
    while iteration < max_iterations
        iteration = iteration + 1;
        found_next = false;
        
        % Look for next boundary pixel in 8-neighborhood
        for d = 1:size(directions, 1)
            next_pos = current + directions(d,:);
            
            % Check if this position is a boundary pixel
            idx = find(rows == next_pos(1) & cols == next_pos(2) & ~visited, 1);
            
            if ~isempty(idx)
                current = next_pos;
                ordered_boundary = [ordered_boundary; current];
                visited(idx) = true;
                found_next = true;
                break;
            end
        end
        
        % If back to start or no more neighbors, stop
        if ~found_next || all(visited)
            break;
        end
    end
    
    % ---- STEP 4: Calculate perimeter ----
    perimeter = 0;
    num_points = size(ordered_boundary, 1);
    
    for k = 1:num_points
        next_k = mod(k, num_points) + 1;
        dx = ordered_boundary(next_k, 2) - ordered_boundary(k, 2);
        dy = ordered_boundary(next_k, 1) - ordered_boundary(k, 1);
        perimeter = perimeter + sqrt(dx^2 + dy^2);
    end
    
    disp(['Perimeter: ', num2str(perimeter)]);
    
    % ---- STEP 5: Calculate circularity ----
    M00 = sum(binaryMask(:));
    circularity = (4*pi*M00)/(perimeter^2);
    disp(['Circularity: ', num2str(circularity)]);
    
    % ---- STEP 6: Visualize ----
    figure, imshow(binaryMask);
    hold on;
    title(['Perimeter: ' num2str(perimeter, '%.2f') ', Circularity: ' num2str(circularity, '%.3f')]);
    hold off;
end

function draw_box_and_centroids(mask)
    if ~any(mask(:)), return; end
    [rows, cols] = find(mask);
    xmin = min(cols); xmax = max(cols);
    ymin = min(rows); ymax = max(rows);

    % bounding box midpoint centroid
    cx_box = (xmin + xmax) / 2;
    cy_box = (ymin + ymax) / 2;

    % moment centroid (note x corresponds to cols, y to rows)
    [X, Y] = meshgrid(1:size(mask,2), 1:size(mask,1));
    mask = double(mask); % ensure mask is numeric 0/1

    M00 = sum(mask(:));
    if M00 == 0
        return;
    end
    M10 = sum(sum(X .* mask));
    M01 = sum(sum(Y .* mask));
    cx_mom = M10 / M00;
    cy_mom = M01 / M00;

    % draw rectangle
    rectangle('Position', [xmin, ymin, xmax - xmin + 1, ymax - ymin + 1], 'EdgeColor', 'y', 'LineWidth', 2);

    % plot centroids
    plot(cx_box, cy_box, 'o', 'Color', 'w', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
    plot(cx_mom, cy_mom, 'x', 'Color', 'r', 'MarkerSize', 10, 'LineWidth', 1);

end

function calculate_orientation(mask)
    [X, Y] = meshgrid(1:size(mask,2), 1:size(mask,1));
    mask = double(mask); % ensure mask is numeric 0/1
    M00 = sum(mask(:));
    if M00 == 0
        return;
    end
    M10 = sum(sum(X .* mask));
    M01 = sum(sum(Y .* mask));
    cx_mom = M10 / M00;
    cy_mom = M01 / M00;
    M20 = sum(sum((X.^2) .* mask));
    M02 = sum(sum((Y.^2) .* mask));
    M11 = sum(sum(X .* Y .* mask));

    u00=M00;
    u01=0;
    u10=0;
    u11 = M11 - cx_mom*M01 - cy_mom*M10 + cx_mom*cy_mom*M00;
    u20 = M20 - 2*cx_mom*M10 + cx_mom^2*M00;
    u02 = M02 - 2*cy_mom*M01 + cy_mom^2*M00;

    J=[u20 u11; u11 u02]
    e=eig(J)

    lambda = sort(e, 'descend');

    
    lambda1 = lambda(1);
    lambda2 = lambda(2);
    
    a=2*sqrt(lambda1/M00);
    b=2*sqrt(lambda2/M00);
    v1 = null(J - lambda1 * eye(2));
    x1 = v1(1);   % first component (x)
    y1 = v1(2);   % second component (y)
    % Calculate the orientation of the bounding box
    theta = atan2(y1, x1);

        t = linspace(0, 2*pi, 200); % parametric angle
    ellipse_x = a * cos(t);
    ellipse_y = b * sin(t);

    % Rotate the ellipse by theta
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rotated_coords = R * [ellipse_x; ellipse_y];

    % Shift ellipse to the centroid (cx_mom, cy_mom)
    ellipse_x_rot = rotated_coords(1, :) + cx_mom;
    ellipse_y_rot = rotated_coords(2, :) + cy_mom;

    % Plot ellipse
    plot(ellipse_x_rot, ellipse_y_rot, 'g', 'LineWidth', 2);

    % Major axis (along θ)
    major_x = [cx_mom - a*cos(theta), cx_mom + a*cos(theta)];
    major_y = [cy_mom - a*sin(theta), cy_mom + a*sin(theta)];
    plot(major_x, major_y, 'r-', 'LineWidth', 2); % red = major axis

    % Minor axis (perpendicular to θ)
    minor_x = [cx_mom - b*sin(theta), cx_mom + b*sin(theta)];
    minor_y = [cy_mom + b*cos(theta), cy_mom - b*cos(theta)];
    plot(minor_x, minor_y, 'b--', 'LineWidth', 2); % blue dashed = minor axis
   
end

    [m, n] = size(binaryMask);
    
    % ---- STEP 1: Find ALL boundary pixels manually ----
    boundary = zeros(m, n);
    
    for i = 2:m-1
        for j = 2:n-1
            if binaryMask(i,j) == 1
                % Check 8-neighbors: if any neighbor is 0, this is a boundary pixel
                neighbors = [binaryMask(i-1,j-1), binaryMask(i-1,j), binaryMask(i-1,j+1), ...
                            binaryMask(i,j-1),                        binaryMask(i,j+1), ...
                            binaryMask(i+1,j-1), binaryMask(i+1,j), binaryMask(i+1,j+1)];
                
                if any(neighbors == 0)
                    boundary(i,j) = 1;
                end
            end
        end
    end
