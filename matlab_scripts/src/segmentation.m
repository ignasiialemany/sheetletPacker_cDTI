function [im_cells, im_groups, im_histology, resolution] = segmentation(tif)
    % Extract the segmentation from a TIF image.
    %   (INPUT)
    %   tif : TIF image with 3 pages (see below)
    %   [OUTPUT]
    %   im_cells     : image of cells (TIF Page 3 [combined with p2])
    %   im_groups    : image of groups (TIF Page 2)
    %   im_histology : original histology image (TIF Page 1)

    %TODO: alternatively, tif can be struct:
    % struct('page1', [], 'page2', [], 'page3', [], 'XResolution', [], 'YResolution', [])

    % meta data
    info = imfinfo(tif);
    info_main = info(1);
    resolution = [info_main.XResolution, info_main.YResolution];
    % add this tag with `exiftool -XResolution=1 -YResolution=2 file.tif`

    % read from TIF
    im_histology = imread(tif, 1); % original histology image
    im_cells = imread(tif, 2); % BW in RGB
    im_groups = imread(tif, 3); % different colours denote different groups

    % remove alpha
    if size(im_histology, 3) == 4
        im_histology = im_histology(:, :, 1:3);
    end
    if size(im_cells, 3) == 4
        im_cells = im_cells(:, :, 1:3);
    end
    if size(im_groups, 3) == 4
        im_groups = im_groups(:, :, 1:3);
    end

    bw_myocytes = imbinarize(rgb2gray(im_cells)); % is already bw
    im_cells = im_groups;
    im_cells(repmat(~bw_myocytes, [1, 1, 3])) = 0;

end
