function main
    % Detects different fruits (apple, capsicum, lemon) from an image using
    % RGB thresholding, extracts the largest component of each fruit, and 
    % visualizes their bounding boxes and centroids.
    
    clc; close all;

    %Read the image
    img = imread('Fruits.jpg');
    %Apply the gaussain filter for smoothing
    img = imgaussfilt(img, 2);

    % Separate RGB channels for thresholding
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);

    % Color Thresholding for each fruit
    % Creating Binary Masks for each fruits using RGB Channels
    apple_mask = (R > 150 & R < 255) & (G > 20 & G < 110) & (B > 60 & B < 120);     
    capsicum_mask = ((R > 150 & R < 255) & (G > 70 & G < 150) & (B < 5));
    lemon_mask = (R > 190 & R < 255) & (G > 150 & G < 255) & (B > 10 & B<120);
    
    % Each mask is converted to double and filters are applied to make the image smooth and remove noise 
    
    apple_mask    = refine_mask(apple_mask);
    figure; imshow(apple_mask); title('Apple Mask');
    
    capsicum_mask = refine_mask(capsicum_mask);
    figure; imshow(capsicum_mask); title('Capsicum Mask');

    lemon_mask    = refine_mask(lemon_mask);
    figure; imshow(lemon_mask); title('Lemon Mask');

    % Extract largest connected component for each fruit
    apple_clean = get_largest_component(apple_mask);
    capsicum_clean = get_largest_component(capsicum_mask);
    lemon_clean = get_largest_component(lemon_mask);

    % Display the results
    figure; imshow(img); hold on;
    draw_box_and_centroids(apple_clean, 'r', 'Apple');
    draw_box_and_centroids(capsicum_clean, 'c', 'Capsicum');
    draw_box_and_centroids(lemon_clean, 'y', 'Lemon');
    hold off;
end


%% ========================================================================
%%                       HELPER FUNCTIONS
%% ========================================================================

function mask_out=refine_mask(mask)
    % Converts a logical mask to double, applies Gaussian and median filtering
    % for smoother edges and less noise.
    % Input:
    %   mask_in  - Binary mask of a segmented fruit
    % Output:
    %   mask_out - Cleaned binary mask after filtering
    mask = double(mask);
    mask = imgaussfilt(mask, 2);        % Smoothen edges
    mask = mask > 0.5;                  % Re-binarize
    mask = medfilt2(mask, [5 5]);       % Median filter to remove salt-and-pepper noise
    mask_out = mask;
end

function largest_mask = get_largest_component(bin_mask)
    % get_largest_component
    % Performs row-by-row connected component labeling (4-connectivity)
    % and returns a logical mask containing only the largest component.
    %
    % Inputs:
    %   bin_mask : logical 2D matrix
    % Output:
    %   largest_mask : logical 2D mask of the largest connected region

    bin_mask = logical(bin_mask);
    [rows, cols] = size(bin_mask);

    % Initialize label matrix
    labels = zeros(rows, cols, 'uint32');
    current_label = 0;

    % Union-Find setup
    parent = uint32(1:rows * cols);  % parent array (upper bound)
    
    % --- Pass 1: assign provisional labels and record equivalences ---
    for r = 1:rows
        for c = 1:cols
            if ~bin_mask(r, c)
                continue;  % background
            end

            % Check left and top neighbors (4-connectivity)
            left_label = 0;
            top_label = 0;

            if c > 1 && labels(r, c-1) > 0
                left_label = labels(r, c-1);
            end
            if r > 1 && labels(r-1, c) > 0
                top_label = labels(r-1, c);
            end

            if left_label == 0 && top_label == 0
                % New component
                current_label = current_label + 1;
                labels(r, c) = current_label;
                parent(current_label) = current_label;
            elseif left_label ~= 0 && top_label == 0
                labels(r, c) = left_label;
            elseif left_label == 0 && top_label ~= 0
                labels(r, c) = top_label;
            else
                % Both neighbors have labels → possible equivalence
                min_label = min(left_label, top_label);
                max_label = max(left_label, top_label);
                labels(r, c) = min_label;
                % Record equivalence
                parent = union_labels(parent, min_label, max_label);
            end
        end
    end

    % --- Pass 2: resolve label equivalences ---
    for r = 1:rows
        for c = 1:cols
            if labels(r, c) > 0
                labels(r, c) = find_root(parent, labels(r, c));
            end
        end
    end

    % --- Compute component areas ---
    unique_labels = unique(labels(labels > 0));
    if isempty(unique_labels)
        largest_mask = false(size(bin_mask));
        return;
    end
    areas = zeros(size(unique_labels), 'uint32');
    for i = 1:numel(unique_labels)
        areas(i) = sum(labels(:) == unique_labels(i));
    end

    % --- Find the largest component ---
    [~, idx] = max(areas);
    best_label = unique_labels(idx);

    % --- Create mask for the largest component ---
    largest_mask = (labels == best_label);
end

%% ---------- helper: union-find "find root" ----------
function root = find_root(parent, x)
    while parent(x) ~= x
        x = parent(x);
    end
    root = x;
end

%% ---------- helper: union-find "union" ----------
function parent = union_labels(parent, a, b)
    ra = find_root(parent, a);
    rb = find_root(parent, b);
    if ra ~= rb
        parent(rb) = ra;
    end
end


%% ---------- helper: draw box + centroids ----------
function draw_box_and_centroids(mask, color, labelStr)
    if ~any(mask(:)), return; end
    [rows, cols] = find(mask);
    xmin = min(cols); xmax = max(cols);
    ymin = min(rows); ymax = max(rows);

    % bounding box midpoint centroid
    cx_box = (xmin + xmax) / 2;
    cy_box = (ymin + ymax) / 2;

    % moment centroid (note x corresponds to cols, y to rows)
    [X, Y] = meshgrid(1:size(mask,2), 1:size(mask,1));
    M00 = sum(mask(:));
    M10 = sum(sum(X .* mask));
    M01 = sum(sum(Y .* mask));
    cx_mom = M10 / M00;
    cy_mom = M01 / M00;

    % draw rectangle
    rectangle('Position', [xmin, ymin, xmax - xmin + 1, ymax - ymin + 1], 'EdgeColor', color, 'LineWidth', 2);

    % plot centroids
    plot(cx_box, cy_box, 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 8);
    plot(cx_mom, cy_mom, 'x', 'Color', color, 'MarkerSize', 10, 'LineWidth', 1);

    % label
    text(cx_box + 5, cy_box - 10, labelStr, 'Color', color, 'FontSize', 11, 'FontWeight', 'bold');
end
