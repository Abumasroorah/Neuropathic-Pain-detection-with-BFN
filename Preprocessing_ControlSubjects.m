% This files contains the code for preprocessing the data for proper wPLI
% computation (for the normal control subjects, the data preprocessing is a crucial stage that also requires
% preprocessing experience as wrongly preprocessed data may go in a long
% way to affect the entire analysis. The code here only considered the eyes
% closed condition and processing the eyes open part will also follow the
% same procedure.

% Note that if the EEG data fails to download automatically, the data could
% be easily downloaded to the local drive and proccessed accordingly.

clear; clc;
eeglab; close;

% Create a directory to store the dataset
dataDir = fullfile(pwd, 'MDD_EEG_Data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

% Define the URL of the dataset
datasetURL = 'https://figshare.com/ndownloader/files/5612735'; % Direct link to the dataset ZIP file

% Define the path to save the downloaded ZIP file
zipFilePath = fullfile(dataDir, 'MDD_EEG_Dataset.zip');

% Download the dataset
websave(zipFilePath, datasetURL);

% Unzip the dataset
unzip(zipFilePath, dataDir);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the control EC data
% List all .mat files in the dataset directory
matFiles = dir(fullfile(dataDir, '*.mat'));

% Initialize a list to store filenames of normal control subjects
normalControlFiles = {};

% Loop through each file to identify normal control subjects with eyes closed data
for i = 1:length(matFiles)
    fileName = matFiles(i).name;
    filePath = fullfile(dataDir, fileName);
    
    % Load the .mat file
    data = load(filePath);
    
    % Check if the subject is a normal control and the data is from eyes closed session
    % This check depends on how the metadata is stored in the .mat files
    % For example, if there's a variable 'group' indicating 'control' and 'condition' indicating 'eyes_closed'
    if isfield(data, 'group') && isfield(data, 'condition')
        if strcmpi(data.group, 'control') && strcmpi(data.condition, 'eyes_closed')
            normalControlFiles{end+1} = fileName; %#ok<SAGROW>
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocess
% Define the duration to extract (e.g., 10 seconds)
segment_duration_sec = 10;

% Loop through each normal control subject
for i = 1:length(normalControlFiles)
    fprintf('Processing Subject %02d...\n', i);
    fileName = normalControlFiles{i};
    filePath = fullfile(dataDir, fileName);
    
    % Load the EEG data
    S = load(filePath);
    EEG = S.EEG; % Adjust this line based on the actual variable name in the .mat file
    
    EEG = eeg_checkset(EEG);
    
    % === Extract Eyes Closed Segment ===
    % Assumption: First event corresponds to Eyes Closed
    start_sample = EEG.event(1).latency;
    end_sample = start_sample + segment_duration_sec * EEG.srate - 1;
    
    % Extract the segment
    EEG = pop_select(EEG, 'point', [start_sample end_sample]);
    EEG = eeg_checkset(EEG);
    
    % === Preprocessing Steps ===
    
    % Notch filter at 60 Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1);
    
    % Bandpass filter between 1–45 Hz
    EEG = pop_eegfiltnew(EEG, 1, 45);
    
    % Resample to 256 Hz
    EEG = pop_resample(EEG, 256);
    
    % Rereference to average
    EEG = pop_reref(EEG, []);
    
    % Line noise removal (requires 'cleanLineNoise' function)
    signal = struct('data', EEG.data, 'srate', EEG.srate);
    lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
        'lineNoiseChannels', 1:EEG.nbchan, ...
        'Fs', EEG.srate, ...
        'lineFrequencies', [60 120 180 240], ...
        'p', 0.01, ...
        'fScanBandWidth', 2, ...
        'taperBandWidth', 2, ...
        'taperWindowSize', 4, ...
        'taperWindowStep', 1, ...
        'tau', 100, ...
        'pad', 2, ...
        'fPassBand', [0 EEG.srate/2], ...
        'maximumIterations', 10);
    [clnOutput, ~] = cleanLineNoise(signal, lineNoiseIn);
    EEG.data = clnOutput.data;
    
    % Replace non-finite values with zeros
    EEG.data(~isfinite(EEG.data)) = 0;
    
    % Clean raw data (requires 'clean_rawdata' plugin)
    originalEEG = EEG;
    EEG = clean_rawdata(originalEEG, 5, -1, 0.85, 4, 20, 0.25);
    
    % Interpolate channels
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
    
    % Run ICA
    EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');
    
    % Remove selected ICA components (manual or automated selection)
    % For example, remove components 2, 5, and 12
    EEG = pop_subcomp(EEG, [2 5 12]);
    
    % Remove baseline
    EEG = pop_rmbase(EEG, []);
    
    % Optional: High-pass filter at 30 Hz
    EEG = pop_eegfiltnew(EEG, [], 30);
    
    % Extract and reshape data
    y = EEG.data;
    y = y([1:19 22:24], 251:end);  % Adjust channels and start point as needed
    outData = double(reshape(y, 22, [], 1));
    
    % Save preprocessed data
    varName = sprintf('out_%02d', i);
    temp.(varName) = outData;
    saveName = sprintf('Sub_%02d_EyesClosed.mat', i + 36); % Labels from Sub_37 to Sub_51
    savePath = fullfile(dataDir, saveName);
    save(savePath, '-struct', 'temp');
    clear temp;
end

fprintf('All subjects'' Eyes Closed data processed and saved.\n');

