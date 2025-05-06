function [bw_stack] = split(bw)
    % create a 3D stack of bw images, one layer per object
    % Does not support GPUArray

    % isolate all 4-connected (avoids joining neighbouring objects) cells
    cc = bwconncomp(bw, 4); % bwconncomp does not work with GPU array

    % size of the problem
    N = cc.NumObjects;
    size_bw = size(bw);

    % loop
    bw_stack = cell(N,1);
    for i = 1:N
        bw_i = zeros(size_bw); % false initially
        bw_i(cc.PixelIdxList{i}) = true; % set pixels of this object as true
        bw_stack{i} = bw_i;
    end

    %{
    for iSheetlet = 1:cc.NumObjects
        % find pixel of one myocyte inside this sheetlet
        [row, col] = find(bw_sheetlet_thick & bw_myocytes, 1); % find first matching pixel
        
        % assign colours
        for c = 1:3
            im_sheetlets_C = im(:, :, c); % read current
            im_sheetlets_C(bw_sheetlet_thick) = im_segmentation(row, col, c); % modify
            im(:, :, c) = im_sheetlets_C; % write modified
        end
        
    end
    %}

end
