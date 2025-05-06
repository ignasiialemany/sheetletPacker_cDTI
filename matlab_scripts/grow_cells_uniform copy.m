function grow_cells_uniform(cells, sheetlet, resx, resy, targetECV, filename)

    data_cells = load(cells);
    bw_cells_resized = data_cells.cells_image;
    data_sheetlet = load(sheetlet);
    bw_groups_resized = data_sheetlet.sheetlet_image;

    addpath("src/");
  
    % Calculate the resolution for the resized images
    resolutionX = resx;
    resolutionY = resy;
    
    for k=1:15
        bw_cells_resized = imerode(bw_cells_resized,strel("disk",1));
    end
    bw_cells_resized = imdilate(bw_cells_resized,strel("disk",10));
    

    %[bw_cells_resized, bw_groups_resized, new_res] = thicken(bw_cells_resized, bw_groups_resized, 1);

    % Define Parameters and Call the Thicken Function
    useGPU = true;
    boundaries_sheetlet = bwboundaries(bw_groups_resized);
    x_sheetlet = boundaries_sheetlet{1}(:, 2) * resolutionX;
    y_sheetlet = (size(bw_cells_resized, 1) - boundaries_sheetlet{1}(:, 1)) * resolutionY;
    centroid_x = (max(x_sheetlet) - min(x_sheetlet)) * 0.5 + min(x_sheetlet);
    centroid_y = (max(y_sheetlet) - min(y_sheetlet)) * 0.5 + min(y_sheetlet);

    % Label cells and remove random cells
    [labeled, totalCells] = bwlabel(bw_cells_resized);

    n_poly_approx = int32((sum(bw_groups_resized(:)) * resolutionY * resolutionX * 1000 * 1000) / 150) - 1;
    n_to_remove = totalCells - n_poly_approx;
    disp(n_to_remove);

    if n_to_remove > 0
        img1 = removeRandomCells(bw_cells_resized, n_to_remove);
        img1 = logical(img1);
        %img1 = bw_cells_resized;
    else
        img1 = bw_cells_resized;
    end

    [labeled, totalCells] = bwlabel(img1);


    % Thicken and group cells
    ecv = 1 - sum(bw_cells_resized(:)) / sum(bw_groups_resized(:));

    % Calculate cell areas
    [labeledCells, numCells] = bwlabel(img1);
    boundaries = bwboundaries(img1);
    cellStats = regionprops(labeledCells, 'Area');
    areas = [];
    for cellIdx = 1:length(cellStats)
        cellArea = cellStats(cellIdx).Area * resolutionX * resolutionY;
        areas = [areas, cellArea];
    end

    % Extract and save polygons
    polygons = [];
    total_area = 0;
    all_areas = [];
    tolerance = min(resolutionX, resolutionY) * 0.5;

    for i = 1:length(boundaries)
        boundary = boundaries{i};
        y = (size(bw_cells_resized, 1) - boundary(:, 1)) * resolutionY;
        x = boundary(:, 2) * resolutionX;
        x_2 = x - centroid_x;
        y_2 = y - centroid_y;
        poly = polyshape(x_2, y_2);

        if poly.area > 10 / (1000 * 1000)
            vertices = poly.Vertices;
            nanRows = any(isnan(vertices), 2);
            vertices = vertices(~nanRows, :);
            if ~isequal(vertices(1, :), vertices(end, :))
                vertices(end + 1, :) = vertices(1, :);
            end
            polygon.vertices = vertices;
            polygons = [polygons, polygon];
            total_area = total_area + poly.area
        end

    end

    filepath = "../matlab_outputs/beforethicken_" + filename + ".mat";
    display(filepath);
    save(filepath,"polygons");

    [img1, img2, new_res] = thicken(img1, bw_groups_resized, 1);


    while ecv >= targetECV
        [img1, img2, new_res] = thicken(img1, img2, 1);
        ecv = 1 - sum(img1(:)) / sum(img2(:));
    end

    if ecv >= targetECV
        disp("We are good");
    end

    % Check for bounding box overlap
    bw_cells_thickened = img1;
    [labeledCells, numCells] = bwlabel(bw_cells_thickened);
    cellStats = regionprops(labeledCells, 'BoundingBox');
    potentialOverlap = false(numCells, numCells);

    for i = 1:numCells
        for j = i + 1:numCells
            bbox1 = cellStats(i).BoundingBox;
            bbox2 = cellStats(j).BoundingBox;
            if bbox1(1) < bbox2(1) + bbox2(3) && bbox1(1) + bbox1(3) > bbox2(1) && ...
                    bbox1(2) < bbox2(2) + bbox2(4) && bbox1(2) + bbox1(4) > bbox2(2)
                potentialOverlap(i, j) = true;
                potentialOverlap(j, i) = true;
            end
        end
    end

    % Check for actual polygon overlap and merge regions
    boundaries = bwboundaries(bw_cells_thickened);
    for i = 1:numCells
        overlappingCells = find(potentialOverlap(i, :));
        for j = overlappingCells
            if i ~= j && potentialOverlap(i, j)
                poly1 = boundaries{i};
                poly2 = boundaries{j};
                mask1 = poly2mask(poly1(:, 2), poly1(:, 1), size(bw_cells_thickened, 1), size(bw_cells_thickened, 2));
                mask2 = poly2mask(poly2(:, 2), poly2(:, 1), size(bw_cells_thickened, 1), size(bw_cells_thickened, 2));
                if any(any(mask1 & mask2))
                    labeledCells(labeledCells == j) = i;
                end
            end
        end
    end

    % Relabel the merged cells and convert back to binary
    mergedCells = bwlabel(labeledCells > 0);
    bw_cells_merged = mergedCells > 0;

    % Plot the results
    figure;
    imshow(bw_cells_merged);
    hold on;
    boundaries = bwboundaries(bw_cells_merged);
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
    end
    hold off;

    % Calculate cell areas
    %[labeledCells, numCells] = bwlabel(bw_cells_thickened);
    %boundaries = bwboundaries(bw_cells_thickened);
    cellStats = regionprops(labeledCells, 'Area');
    areas = [];
    for cellIdx = 1:length(cellStats)
        cellArea = cellStats(cellIdx).Area * resolutionX * resolutionY;
        areas = [areas, cellArea];
    end

    % Extract and save polygons
    polygons = [];
    total_area = 0;
    all_areas = [];
    tolerance = min(resolutionX, resolutionY) * 0.5;

    for i = 1:length(boundaries)
        boundary = boundaries{i};
        y = (size(bw_cells_resized, 1) - boundary(:, 1)) * resolutionY;
        x = boundary(:, 2) * resolutionX;
        x_2 = x - centroid_x;
        y_2 = y - centroid_y;
        poly = polyshape(x_2, y_2);

        if poly.area > 10 / (1000 * 1000)
            vertices = poly.Vertices;
            nanRows = any(isnan(vertices), 2);
            vertices = vertices(~nanRows, :);
            if ~isequal(vertices(1, :), vertices(end, :))
                vertices(end + 1, :) = vertices(1, :);
            end
            polygon.vertices = vertices;
            polygons = [polygons, polygon];
            total_area = total_area + poly.area
        end

    end
    filepath = "../matlab_outputs/" + filename + ".mat";
    display(filepath);
    save(filepath,"polygons");
end


    