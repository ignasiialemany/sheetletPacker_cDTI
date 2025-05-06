function [] = write_poly(filename, polygons)

    fid = fopen(filename, 'w');

    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'write_poly.m\n'); %TODO: make optional
    fprintf(fid, 'ASCII\n'); %TODO: make optional

    polygons = polygons(:); % code below requires this in column format

    poly_vert = arrayfun(@(p) p.Vertices, polygons, 'UniformOutput', false);
    poly_size = cellfun(@(v) size(v, 1), poly_vert);
    vertices_xy = cell2mat(poly_vert); % concatenate

    N_poly = numel(poly_vert); % same as poly_vert
    N_vert = size(vertices_xy, 1);

    vertices = [vertices_xy, zeros(N_vert, 1)]; % pad with z=0

    fprintf(fid, 'DATASET POLYDATA\n');
    fprintf(fid, 'POINTS %d %s\n', N_vert, 'double');
    fprintf(fid, '%g %g %g\n', transpose(vertices));
    fprintf(fid, 'POLYGONS %d %d\n', N_poly, sum(1+poly_size));
    cumsum_poly_size = cumsum([0; poly_size]);
    for i = 1:N_poly
        Ni = poly_size(i);
        fmt = ['%d', repmat(' %d', [1, Ni]), '\n'];
        connectivities = cumsum_poly_size(i) + (0:Ni-1); % 0-based indexing
        fprintf(fid, fmt, Ni, connectivities);
    end

    fclose(fid);

end
