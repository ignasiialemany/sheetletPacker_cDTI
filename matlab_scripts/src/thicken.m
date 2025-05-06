function [bw_cells, bw_groups, nIter] = thicken(bw_cells, bw_groups, nIterTarget, useGPU)
    %THICKEN Thicken cells N times or until no further change occurs
    %   (INPUT)
    %   bw_cells  : binary image, representing the cells
    %   bw_groups : binary image, representing the groups of cells
    %   nIter     : number of iterations
    %   [OUTPUT]
    %   bw_cells  : binary image, representing the morphed cells
    %   bw_groups : binary image, representing the morphed groups of cells
    %   nIter     : number of iterations (may be less than input)

    if nIterTarget == 0
        nIter = 0;
        return % no need for silly overheads
    end

    % Use of GPU considerably accelerates this
    if nargin < 4
        useGPU = true; % use GPU by default, but can be overruled
    end
    bw_cells = gpu(bw_cells);
    bw_groups = gpu(bw_groups);

    %TODO: groups should be able to morph too

    % grow objects
    nIter = 0;

    outputFolder = 'thickening_iterations';
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    while nIter < nIterTarget

        % store previous version
        bw_morphed_old = bw_cells;

        % morph (grow)
        bw_cells = bwmorph(bw_cells, 'spur', inf()); % remove spurious pixels
        bw_cells = bwmorph(bw_cells, 'thicken', 1);

        % clip to group
        bw_cells = bw_cells & bw_groups;
        bw_cells = bwmorph(bw_cells, 'spur', inf()); % remove spurious pixels

        % test for change
        nochange = isequal(bw_cells, bw_morphed_old);
        if nochange
            nIter = 0;
            break % end looping
        else
            nIter = nIter + 1; % increment
        end

         % Plotting the current state of bw_cells
        %figure('visible', 'off');
        %imshow(bw_cells);
        %title(sprintf('Iteration %d', nIter));
        %saveas(gcf, fullfile(outputFolder, sprintf('Iteration_%d.png', nIter)));
        %close;
    end

    % Return as CPU data
    bw_cells = cpu(bw_cells);
    bw_groups = cpu(bw_groups);


    %%% Auxiliary functions for CPU/GPU data management

    function [dat] = cpu(dat)
        % Send data to the CPU.
        if has_GPU()
            dat = gather(dat);
        end
    end

    function [dat] = gpu(dat)
        % Send data to the GPU.
        if has_GPU()
            dat = gpuArray(dat);
        end
    end

    function [tf] = has_GPU()
        % Test if GPU is available.
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