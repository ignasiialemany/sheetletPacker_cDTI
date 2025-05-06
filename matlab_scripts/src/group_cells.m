function [bw_sheetlet] = thicken_sheetlet(bw_sheetlet, bw_sheetlet_total_area, minGroupSize, resx, resy)

% Label connected components
[labeledCells, numCells] = bwlabel(bw_sheetlet);
numCellsOriginal = numCells;

% Get properties of labeled regions
cellProps = regionprops(labeledCells, 'BoundingBox', 'PixelIdxList');

% Initialize cell array to hold group images and their masks
groupsWithMasks = {};

usedCells = false(numCells, 1);

% Group cells based on neighboring proximity
groupLabel = zeros(size(bw_sheetlet));  % Initialize the label matrix for groups
currentGroupNum = 1;  % Initialize group number

% Define the threshold for bounding box overlap ratio
overlapThreshold = 0.08;  % Adjust this threshold as needed
boundingBoxMargin = 5;  % Adjust the margin as needed

% First pass to process groups of cells
for i = 1:numCells
    if usedCells(i)
        continue;
    end
    
    bw_check_areas = imerode(bw_sheetlet,strel("disk",4));
    [labeledCells, numCells] = bwlabel(bw_check_areas);
    if numCellsOriginal ~= numCells
        x=3;
    end
    

    % Get properties of labeled regions
    cellProps = regionprops(labeledCells, 'BoundingBox', 'PixelIdxList');


    % Initialize group with the current cell
    currentGroup = i;
    usedCells(i) = true;
    
    % Get the current cell's bounding box and expand it
    bbox1 = cellProps(i).BoundingBox;
    expandedBbox1 = [bbox1(1) - boundingBoxMargin, bbox1(2) - boundingBoxMargin, bbox1(3) + 2*boundingBoxMargin, bbox1(4) + 2*boundingBoxMargin];
    
    % Get the pixel list of the current cell
    pixelList = cellProps(i).PixelIdxList;

    % Check the neighbors of the current cell
    for j = 1:numCells
        if i == j || usedCells(j)
            continue;
        end
        
        % Get the bounding box of the neighboring cell and expand it
        bbox2 = cellProps(j).BoundingBox;
        expandedBbox2 = [bbox2(1) - boundingBoxMargin, bbox2(2) - boundingBoxMargin, bbox2(3) + 2*boundingBoxMargin, bbox2(4) + 2*boundingBoxMargin];
        
        % Calculate the overlap ratio using the expanded bounding boxes
        overlapRatio = bboxOverlapRatio(expandedBbox1, expandedBbox2, 'union');
        
        % Check if the overlap ratio meets the threshold
        if overlapRatio > 0
            currentGroup = [currentGroup; j];
            usedCells(j) = true;
        end
    end

    % Check if the group meets the minimum group size requirement
    if numel(currentGroup) >= minGroupSize && numel(currentGroup) < 5
        % Create a new logical image for the group
        groupImage = ismember(labeledCells, currentGroup);
        imageSize = size(groupImage);

        % Determine the convex hull of the group
        boundaries = bwboundaries(groupImage);
        allPoints = vertcat(boundaries{:});
        k = convhull(allPoints(:,2), allPoints(:,1));
        convexHullPoints = allPoints(k, :);

        % Create a binary mask for the convex hull
        hullPolyshape = polyshape(convexHullPoints(:,2), convexHullPoints(:,1));
        [columnsInImage, rowsInImage] = meshgrid(1:size(groupImage, 2), 1:size(groupImage, 1));
        inHull = inpolygon(columnsInImage, rowsInImage, hullPolyshape.Vertices(:,1), hullPolyshape.Vertices(:,2));
        binaryMask = false(size(groupImage));
        binaryMask(inHull) = true;
        
        % Iterate over each pixel in the groupImage
        for row = 1:imageSize(1)
            for col = 1:imageSize(2)
                % Check if the current pixel is inside the convex hull
                if binaryMask(row, col) == 1 && ~ismember(labeledCells(row, col), currentGroup) && labeledCells(row, col) ~= 0
                    binaryMask(row, col) = 0;
                end
            end
        end

        % Store the group image and the corresponding mask
        binaryMask = imerode(binaryMask, strel("disk", 5));
        groupsWithMasks{end+1} = struct('groupImage', groupImage, 'binaryMask', binaryMask);
        
        close all;
        figure();
        imshowpair(binaryMask,groupImage);

        % Thicken the group to achieve the desired density
         random_packing = rand()*0.05+0.05;
        [cells_thicken, group_thicken, iter] = thicken(groupImage, binaryMask, 1, true);

        area_cells_initial = sum(cells_thicken(:)) * resx * resy;
        total_area = sum(group_thicken(:)) * resx * resy;
        initial_ecs = 1 - (area_cells_initial / total_area);

        while initial_ecs >= random_packing
            [cells_thicken, group_thicken, iter] = thicken(cells_thicken, group_thicken, 1, true);
            area_cells_initial = sum(cells_thicken(:)) * resx * resy;
            total_area = sum(group_thicken(:)) * resx * resy;
            initial_ecs = 1 - (area_cells_initial / total_area);
        end
        
        cells_thicken = imerode(cells_thicken, strel("disk", 1));

        bw_sheetlet = logical(cells_thicken + bw_sheetlet);
        
        close all;
        figure();
        imshow(bw_sheetlet);

        total_ics = sum(bw_sheetlet(:)) * resx * resy * 1000 * 1000;
        total_sheetlet_area = sum(bw_sheetlet_total_area(:)) * resx * resy * 1000 * 1000;
        display(1 - total_ics / total_sheetlet_area);

        % Label the group in the label matrix
        groupLabel(binaryMask) = currentGroupNum;
        currentGroupNum = currentGroupNum + 1;
    end
end

end
