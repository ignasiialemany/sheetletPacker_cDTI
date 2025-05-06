function [im] = colorize(bw_morphed, im_orig)

    bw_orig = im2bin(im_orig);

    % colourise
    [row, col] = find(bw_morphed & bw_orig, 1); % find first common pixel
    colour = im_orig(row, col, :); % colour of this object

    bw3D = int8(repmat(bw_morphed, [1, 1, 3])); % use int8 as temp
    bw3D = reshape(typecast(bw3D(:), class(colour)), size(bw3D));
    im = bw3D .* colour;

end
