function BC = betti_curve(persistence_diagram)
    % Compute Betti curve
    % Assuming persistence_diagram is a cell array of persistence diagrams

    BC = zeros(1, length(persistence_diagram));

    for i = 1:length(persistence_diagram)
        BC(i) = length(persistence_diagram{i}); % Number of points in each persistence diagram
    end
end
