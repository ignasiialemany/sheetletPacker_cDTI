function [new_cells] = removeRandomCells(cells, n_remove)

    % We get labeled image
    [new_cells, numCells] = bwlabel(cells);

    % Get unique cell indices
    index_cells = unique(new_cells);

    % Remove the background index (0)
    index_cells(index_cells == 0) = [];

    % Ensure n_remove does not exceed the number of cells
    n_remove = min(n_remove, length(index_cells));

    % Randomly select unique indices to remove
    indices_to_remove = index_cells(randperm(length(index_cells), n_remove));

    % Create a mask for the cells to be removed
    mask = ismember(new_cells, indices_to_remove);

    % Set the selected cells to 0
    new_cells(mask) = 0;

end