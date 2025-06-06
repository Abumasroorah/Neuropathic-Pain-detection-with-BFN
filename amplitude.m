function AD = amplitude(persistence_diagram)
    % Compute amplitude
    % Assuming persistence_diagram is a cell array of persistence diagrams

    AD = zeros(1, length(persistence_diagram));

    for i = 1:length(persistence_diagram)
        distances = abs(persistence_diagram{i}(:, 2) - persistence_diagram{i}(:, 1));
        AD(i) = norm(distances); % Compute L1 norm
    end
end
