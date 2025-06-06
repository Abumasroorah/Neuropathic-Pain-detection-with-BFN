clear; 
clc;

% Define data directory and sampling frequency
dataDir = fullfile(pwd, 'NeuropathicPainEEG');
Fs = 256;            % Sampling rate
L = 6 * Fs;          % 6 seconds = 1536 samples

% Initialize output
DATA = cell(36,1);

% Loop through all subjects
for j = 1:36
    % Load subject file
    fileName = fullfile(dataDir, sprintf('Sub_%02d.mat', j));
    s = load(fileName);
    
    % Extract variable dynamically (e.g., out_01, out_02, etc.)
    varName = sprintf('out_%02d', j);
    Data = s.(varName);
    
    % Trim data to multiple of 6s and segment
    totalSamples = size(Data, 2);
    numEpochs = floor(totalSamples / L);
    data = Data(:, 1:(numEpochs * L));
    
    % Segment into 6-sec epochs
    AB_data = zeros(numEpochs, size(data,1), L);
    for i = 1:numEpochs
        AB_data(i,:,:) = data(:, (i-1)*L+1 : i*L);
    end
    
    % Store segmented data
    DATA{j} = AB_data;
end

% Save the result
save('PainData_6s.mat', 'DATA2', '-v7.3');
fprintf('Segmentation complete. Data saved as PainData_6s.mat\n');



%  The same process as above should be done for the controls subjects and
%  the NP and Control;s subjects data could be concatenated as DATA2.
