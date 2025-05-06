function [bw_cells, bw_groups, nIter] = thicken_keep(bw_cells, bw_groups, nIterTarget, useGPU, tArea, res_x, res_y)
    %THICKEN Thicken cells N times or until no further change occurs

    % Return immediately if no thickening is requested
    if nIterTarget == 0
        nIter = 0;
        return
    end

    % Use GPU if available
    if nargin < 4
        useGPU = true;
    end
    bw_cells = gpu(bw_cells);
    bw_groups = gpu(bw_groups);

    % Thicken the smaller cells
    nIter = 0;
    while nIter < nIterTarget
        % Separate cells based on area
        CC = bwconncomp(bw_cells);
        stats = regionprops(CC, 'Area');
        bw_large_cells = false(size(bw_cells));
        bw_small_cells = false(size(bw_cells));
        for i = 1:numel(stats)
            if (stats(i).Area)*res_x*res_y > tArea
                bw_large_cells(CC.PixelIdxList{i}) = true;
            else
                bw_small_cells(CC.PixelIdxList{i}) = true;
            end
        end

        %Add padding to large cells
        bw_large_cells_buffer = bwperim(bw_large_cells) + bw_large_cells;
        [r, c] = find(bw_large_cells_buffer);
        valid_indices = r > 1 & r < size(bw_large_cells_buffer, 1) & c > 1 & c < size(bw_large_cells_buffer, 2);
        r = r(valid_indices);
        c = c(valid_indices);
        for k = 1:numel(r)
            bw_large_cells_buffer(r(k)-1:r(k)+1, c(k)-1:c(k)+1) = true;
        end

        % Subtract the buffer zone from the smaller cells
        bw_small_cells = bw_small_cells & ~bw_large_cells;

        bw_morphed_old = bw_small_cells;

        bw_small_cells = bwmorph(bw_small_cells, 'spur', inf());
        bw_small_cells = bwmorph(bw_small_cells, 'thicken', 1);

        bw_small_cells = bw_small_cells & bw_groups & ~bw_large_cells_buffer;
        bw_small_cells = bwmorph(bw_small_cells, 'spur', inf());
        
        bw_cells = bw_large_cells | bw_small_cells;

        if isequal(bw_small_cells, bw_morphed_old)
            break
        else
            nIter = nIter + 1;
        end
    end
    % Combine the cells back together
    
    % Convert back to CPU data
    bw_cells = cpu(bw_cells);
    bw_groups = cpu(bw_groups);

    % Auxiliary functions for CPU/GPU data management

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
                tf = gpuDeviceCount();
            catch
                tf = false;
            end
        end
    end

end
