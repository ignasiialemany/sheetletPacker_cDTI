function [polys] = im2polys(bw)
    % Extract polygons representing the region boundaries in a binary image.

    bw_stack = split(bw);
    polys = cellfun(@(bw) im2poly(bw), bw_stack, 'ErrorHandler', @errfun);

    % remove empty polygons
    polys = polys(arrayfun(@(p) ~isequal(p, polyshape()), polys));

end

function [poly] = errfun(err, varargin)
    % ErrorHandler, so cellfun does not rethrow. We simply want to create an empty polyshape
    poly = polyshape(); % empty one

    %TODO: could recast error as warning?
    %{
    warning(err.message);
    warning(err.identifier, err.message);
    fprintf('error at index %d\n', err.index);
    %}

end

function [poly] = im2poly(bw)
    % Convert a binary image containing a single region to a polygon.

    % get 8-connected (to avoid staircase) boundary coordinates
    boundaries = transpose(bwboundaries(bw, 8, 'noholes'));
    if numel(boundaries) ~= 1
        error('CellMorphing:im2polys:num_boundaries', ...
            'Was expecting one boundary but found %d', numel(boundaries));
    end
    boundary1 = boundaries{1}; % pick the only one

    % convert
    yx_ccw = boundary1; % [y(:),x(:)], anticlockwise
    xy_ccw = fliplr(yx_ccw); % [x(:), y(:)]
    xy = flipud(xy_ccw); %#ok<FLUDLR> {flip better than rot} % clockwise

    % construct
    poly = make_polyshape(xy);
    verify_poly(poly);

end

function [poly] = make_polyshape(xy)
    % Make a polyshape object. Handles exceptions and may throw error.

    msgID1 = 'MATLAB:polyshape:repairedBySimplify'; % intersection or duplicate (collinear) vertices
    msgID2 = 'MATLAB:polyshape:boundary3Points'; % too few points
    warnstruct_old = [warning('query', msgID1), warning('query', msgID2)]; % use "" not ''
    errstate = {'error', 'error'};
    warnstruct_err = warnstruct_old;
    [warnstruct_err.state] = errstate{:};
    warning(warnstruct_err); % https://undocumentedmatlab.com/articles/trapping-warnings-efficiently

    try

        poly = polyshape_(xy); % automatically removes collinear points
        warning(warnstruct_old); % reset upon success

    catch exception % we encountered a warning (or actual error) in @polyshape_

        warning(warnstruct_old); % reset here too

        if strcmp(exception.identifier, msgID1) % ill-defined polygon

            % re-do what caused the warning/error, but be quiet
            warning('off', msgID1); % temporarily turn off
            poly_n = polyshape_(xy);
            warning(warnstruct_old(1)); % reset

            % now diagnose based on constructed polyshape (most likely due to self-intersection)
            if poly_n.NumRegions > 1 % don't use numboundaries, because there might be holes
                poly_n_sorted = poly_n.sortregions('area', 'descend');
                [x, y] = boundary(poly_n_sorted, 1); % extract first (=largest after sorting)
                poly = polyshape_(x, y); % recreate polygon for only this region
            else % something else, but problem was fixed and a single polyshape was created

                % must have been duplicate vertices
                % could also have been that a line was provided
                poly = poly_n; % use already-fixed polyshape
                % either way, if something invalid was generated, it will be caught later

            end

        elseif strcmp(exception.identifier, msgID2) % too few points

            poly = polyshape(); % don't rethrow, just generate empty (will be caught)

        else % unknown error

            rethrow(exception);

        end
    end

end

function [poly] = polyshape_(varargin)
    % Construct polyshape from either (xy) or (x, y)
    poly = polyshape(varargin{:}, 'KeepCollinearPoints', false);
end

function verify_poly(poly)
    % Verify integrity of polygon

    %TODO: do we need to handle inputs with holes? if so, is this even a problem?
    assert(poly.NumHoles==0, ...
        'CellMorphing:im2polys:has_holes', ...
        'Polygon has holes');

    %TODO: should not happen, but test anyway?
    % NaN in vertices only (?) occurs for multiple regions, which we caught.
    assert(~any(isnan(poly.Vertices),'all'), ...
        'CellMorphing:im2polys:nan_vertex', ...
        'Vertices with NaN detected');

    if poly.isequal(polyshape()) % empty
        % likely due to providing a point or a straight line
        error('CellMorphing:im2polys:empty_polygon', ...
            'Could not construct a polygon');
    end

end
