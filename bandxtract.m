function segmentedData = bandxtract(data)
% Assuming your 4D matrix is named 'data' with dimensions 26x27x45x200
% This function xtract the various bands networks from the full band
% Abdul 2023

% Define the frequency bands
freqBands = [1 4; 4 8; 8 13; 13 30; 30 45];

% Initialize a cell array to store segmented data
segmentedData = cell(1, size(freqBands, 1));

% Loop through each frequency band
for i = 1:size(freqBands, 1)
    % Extract frequency range for the current band
    startFreq = freqBands(i, 1);
    endFreq = freqBands(i, 2);
    
    % Find indices corresponding to the frequency range
    indices = startFreq:endFreq;
    
    % Segment the data for the current band
    segmentedData{i} = data(:, :, indices, :);
%     segmentedData{i} = data(:, :, indices, :);
    
end

% Access the segmented data using segmentedData{1}, segmentedData{2}, etc.
% Each cell contains a 4D matrix representing the data for a specific frequency band
