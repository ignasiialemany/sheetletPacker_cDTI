function [bw_cells, bw_groups, nIter] = thicken_area(bw_cells, bw_groups, useGPU, resx, resy, targetECV)

    sheetlet_area_total = sum(bw_groups(:)) * resx * resy * 1000 * 1000;
    if nargin < 3
        useGPU = true;
    end

    if nargin < 6
        error('All inputs (bw_cells, bw_groups, resx, resy, targetECV) must be provided.');
    end

    % Initialize locked areas and original locked cells
    lockedAreas = false(size(bw_cells));
    original_locked_cells = false(size(bw_cells));

    % GPU handling
    bw_cells = gpu(bw_cells);
    bw_groups = gpu(bw_groups);

    outputFolder = 'thickening_iterations';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Get boundaries and create polyshapes
    boundaries = bwboundaries(bw_cells);
    polyshapes = cellfun(@(b) polyshape(b(:,2), b(:,1)), boundaries, 'UniformOutput', false);
    numCells = numel(polyshapes);
    
    % Get the areas of the polyshapes
    startingAreas = cellfun(@area, polyshapes);
    startingAreas = startingAreas * resx * resy * 1000 * 1000;

    indices_fillers = [];
    indices_cells = [];

    for i = 1:length(startingAreas)
        if startingAreas(i) < 25
            indices_fillers = [indices_fillers, i];
        else
            indices_cells = [indices_cells, i];
        end
    end
    cellAreaTargets = getAreasfromPdf(numCells);
    maxAreaAfterThickening = sum(cellAreaTargets);
    available_area = sum(bw_groups(:)) * resx * resy * 1000 * 1000;
    targeted_area = (1 - targetECV) * available_area;

    cellAreaTargets(indices_fillers) = 2;

    ratio_available_targeted = available_area / targeted_area;
    [rows, cols] = size(bw_groups);
    ratio_available_targeted = 1.0;
    while ratio_available_targeted >= 1.05
        subtraction = bw_groups - bw_cells;
        while true
            [rowInds, colInds] = find(subtraction);

            randIndex = randi(length(rowInds));
            centerX = colInds(randIndex);
            centerY = rowInds(randIndex);

            radius = randi([20, 40]);

            [X, Y] = meshgrid(1:cols, 1:rows);

            distanceFromCenter = sqrt((X - centerX).^2 + (Y - centerY).^2);

            circleMask = distanceFromCenter <= radius;
            ones_in_mask = sum(subtraction(circleMask));

            if ones_in_mask < ceil(length(subtraction(circleMask)) * 0.7)
                mask_indices = find(circleMask);
                for idx = mask_indices'
                    if subtraction(idx)
                        bw_groups(idx) = 0;
                    end
                end
                break;
            end
        end

        available_area = sum(bw_groups(:)) * resx * resy * 1000 * 1000;
        ratio_available_targeted = available_area / targeted_area;
    end

    available_area = sum(bw_groups(:)) * resx * resy * 1000 * 1000;
    target_area_ics = (1 - targetECV) * available_area;
    cellAreaTargets  = normrnd(188,45,[numCells,1]);

    cellAreaTargets = cellAreaTargets / (1000 * 1000);
    isAreaLocked = startingAreas < -1;
    ics_area_cells = sum(bw_cells(:)) * resx * resy * 1000 * 1000;
    sheetlet_area = sum(bw_groups(:)) * resx * resy * 1000 * 1000;

    ecv = 1 - ics_area_cells / sheetlet_area;
    nIter = 1;
    oldEcv = ecv;

    while ecv >= targetECV
        
        display(ecv);
        bw_checkareas = imerode(bw_cells, strel("disk", 1));
        
        % Skeletonize the binary image
        % Label the skeletonized image
        [boundaries, labeledCells] = bwboundaries(bw_checkareas);

        % Get the properties of the labeled regions
        polyshapes = cellfun(@(b) polyshape(b(:,2), b(:,1)), boundaries, 'UniformOutput', false);
        centroidsX = zeros(length(boundaries));
        centroidsY = zeros(length(boundaries));

        for j=1:length(polyshapes)
            [centroidsX(j),centroidsY(j)] = centroid(polyshapes{j});
        end
       
        new_locked_area = false(size(bw_cells));

        if nIter > 1
            labelMap = mapCentroids(prevCentroidsX, prevCentroidsY, centroidsX, centroidsY);
        else
            labelMap = 1:numCells;
        end
        prevCentroidsX = centroidsX;
        prevCentroidsY = centroidsY;

        new_areas = [];

        for i = 1:numCells
            mappedIndex = labelMap(i);
            cellArea = area(polyshapes{mappedIndex}) * resx * resy;
            if ismember(mappedIndex, indices_cells)
                new_areas = [new_areas, cellArea];
            end

            if isAreaLocked(mappedIndex) == true
                continue;
            end
            if cellArea >= cellAreaTargets(mappedIndex) && ~isAreaLocked(mappedIndex)
                new_locked_area = new_locked_area | (labeledCells == mappedIndex);
                isAreaLocked(mappedIndex) = true;
            end
        end

        original_locked_cells = new_locked_area | original_locked_cells;
        new_locked_area = bwmorph(new_locked_area, 'thicken', 2);
        lockedAreas = new_locked_area | lockedAreas;

        bw_grow = bw_cells & ~lockedAreas;
        bw_grow = bwmorph(bw_grow, 'spur', inf);
        bw_grow = bwmorph(bw_grow, 'thicken', 1);

        bw_cells = bw_grow & bw_groups & ~lockedAreas;
        bw_cells = bw_cells | original_locked_cells;
        bw_cells = bwmorph(bw_cells, 'spur', inf);

        ics_area_cells = sum(new_areas) * 1000 * 1000;
        ecv = 1 - (ics_area_cells / sheetlet_area_total);
        

        if abs(oldEcv-ecv)<0.0001
            break
        end

        if mod(nIter, 10) == 0
            figure(1);
            imshow(bw_cells);
            title('Cells');

            figure(2)
            cla;
            histogram(new_areas, 50, 'FaceColor', 'b');
            hold on;
            histogram(cellAreaTargets, 50, 'FaceColor', 'r');
            hold off;
            title('Area Histograms');

            drawnow;
        end

        nIter = nIter + 1;
        oldEcv = ecv;
    end

    bw_cells = cpu(bw_cells);
    bw_groups = cpu(bw_groups);

    function [dat] = cpu(dat)
        if has_GPU()
            dat = gather(dat);
        end
    end

    function [dat] = gpu(dat)
        if has_GPU()
            dat = gpuArray(dat);
        end
    end

    function [tf] = has_GPU()
        tf = false;
        if useGPU
            try
                tf = gpuDeviceCount() > 0;
            catch
                tf = false;
            end
        end
    end
end

function labelMap = mapCentroids(prevCentroidsX, prevCentroidsY, currentCentroidsX, currentCentroidsY)
    distances = pdist2([prevCentroidsX(:), prevCentroidsY(:)], [currentCentroidsX(:), currentCentroidsY(:)]);
    [~, labelMap] = min(distances, [], 1);
end

function getAreas = getAreasfromPdf(numCells)
    data = load("areas.mat");
    areas_jan = data.areas_jan;

    bin_edges = 30:0.5:550;

    [counts, edges] = histcounts(areas_jan, bin_edges);

    pdf = counts / sum(counts);

    bin_centers = edges(1:end-1) + diff(edges) / 2;

    cdf = cumsum(pdf);
    [unique_cdf, unique_idx] = unique(cdf);
    unique_bin_centers = bin_centers(unique_idx);

    random_numbers = rand(numCells, 1);

    getAreas = interp1(unique_cdf, unique_bin_centers, random_numbers);
end
