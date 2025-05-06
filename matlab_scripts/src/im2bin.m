function [bw] = im2bin(im)
    % Convert a RGB image into a BW image
    bw = any(im > 0, 3);
end
