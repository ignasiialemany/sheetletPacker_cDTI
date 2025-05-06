function [bw_cells, bw_groups] = shrink(bw_cells, bw_groups, nIter)
    %SHRINK Shrink cells N times
    %   (INPUT)
    %   bw_cells  : binary image, representing the cells
    %   bw_groups : binary image, representing the groups of cells
    %   nIter     : number of iterations
    %   [OUTPUT]
    %   bw_cells  : binary image, representing the morphed cells
    %   bw_groups : binary image, representing the morphed groups of cells

    if nIter == 0
        return % no need for silly overheads
    end

    %TODO: check actual nIter and return it

    bw_cells  = shrink_(bw_cells , nIter(1));
    bw_groups = shrink_(bw_groups, nIter(2));

end

function [bw] = shrink_(bw, nIter)
    % Shrink all objects in an image.

    % we need to operate on each object separately
    bw_stack = split(bw); % @split does not accept GPUArray
    bw_stack = cellfun(@(b) shrink__(b, nIter), bw_stack, 'UniformOutput', false);
    bw = any(cell2mat(reshape(bw_stack, 1,1,[])), 3); % collapse

end

function [bw] = shrink__(bw, nIter)
    % Shrink one object.
    % Does not accept GPUArray

    % create one pixel gap (necessary to properly remove the edge spur)
    bw = padarray(bw, [1, 1], false, 'both');

    % morph
    bw = bwmorph(bw, 'shrink', nIter); % shrink to a POINT
    bw = bwmorph(bw, 'spur', inf()); % remove spurious pixels

    % needed in case the shrunk myocyte is no longer 4-connected
    cc = bwconncomp(bw, 4);
    regions = regionprops(cc, 'Area', 'PixelIdxList');
    [~, id] = max([regions.Area]);
    bw(:) = false;
    bw(regions(id).PixelIdxList) = true;

    % continue
    bw = bwmorph(bw, 'spur', inf()); % remove spurious pixels
    bw = bw(2:end-1, 2:end-1); % remove extra pixel gap

end
