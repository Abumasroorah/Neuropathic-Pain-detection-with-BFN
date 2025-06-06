% This files contains the code for preprocessing the data for proper wPLI
% computation, the data preprocessing is a crucial stage that also requires
% preprocessing experience as wrongly preprocessed data may go in a long
% way to affect the entire analysis. The code here only considered the eyes
% closed condition and processing the eyes open part will also follow the
% same procedure.

% Note that if the EEG data fails to download automatically, the data could
% be easily downloaded to the local drive and proccessed accordingly.

clear; clc;

% Initialize EEGLAB
eeglab; close;

% Create folder to store data
dataDir = fullfile(pwd, 'NeuropathicPainEEG');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

% Step 1: Download and unzip dataset from Mendeley
zipUrl = 'https://data.mendeley.com/public-files/datasets/yj52xrfgtz/4/files/ef963893-2be2-4f35-8132-fd08b3d18640/file_downloaded';
zipFile = fullfile(dataDir, 'neuropathic_pain_data.zip');
websave(zipFile, zipUrl);
unzip(zipFile, dataDir);

% Get all subject .mat files
matFiles = dir(fullfile(dataDir, '*.mat'));
matFiles = natsortfiles({matFiles.name});  % Optional: requires FileExchange function

% Duration to extract in seconds (e.g., 10 seconds)
segment_duration_sec = 10;

% Loop over subjects
for i = 1:length(matFiles)
    fprintf('Processing Subject %02d...\n', i);
    filePath = fullfile(dataDir, matFiles{i});
    
    % Load data
    S = load(filePath);
    EEG = S.EEG;

    EEG = eeg_checkset(EEG);

    % === Extract Eyes Closed Segment ===
    % Assumption: First event is Eyes Closed
    start_sample = EEG.event(1).latency;
    end_sample = start_sample + segment_duration_sec * EEG.srate - 1;

    % Extract only the Eyes Closed segment
    EEG = pop_select(EEG, 'point', [start_sample end_sample]);
    EEG = eeg_checkset(EEG);

    % === Preprocessing ===

    % Notch filter at 60 Hz
    EEG = pop_eegfiltnew(EEG, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1);

    % Bandpass filter 1–45 Hz
    EEG = pop_eegfiltnew(EEG, 1, 45);

    % Resample to 256 Hz
    EEG = pop_resample(EEG, 256);

    % Rereference to average
    EEG = pop_reref(EEG, []);

    % Line noise removal
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

    % Replace non-finite values
    EEG.data(~isfinite(EEG.data)) = 0;

    % Clean raw data
    originalEEG = EEG;
    EEG = clean_rawdata(originalEEG, 5, -1, 0.85, 4, 20, 0.25);

    % Interpolate channels
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

    % Run ICA
    EEG = pop_runica(EEG, 'extended', 1, 'interupt', 'on');

    % Remove selected ICA components (you may automate this)
    EEG = pop_selectcomps(EEG, [1:19 22:24]);
    EEG = pop_subcomp(EEG, [2 5 12]);

    % Remove baseline
    EEG = pop_rmbase(EEG, []);

    % High-pass filter again (optional)
    EEG = pop_eegfiltnew(EEG, [], 30);

    % Extract and reshape
    y = EEG.data;
    y = y([1:19 22:24], 251:end);  % adjust start point if needed
    outData = double(reshape(y, 22, [], 1));

    % Save preprocessed Eyes Closed data
    varName = sprintf('out_%02d', i);
    temp.(varName) = outData;
    saveName = sprintf('Sub_%02d_EyesClosed.mat', i);
    savePath = fullfile(dataDir, saveName);
    save(savePath, '-struct', 'temp');
    clear temp;
end

fprintf('All subjects'' Eyes Closed data processed and saved.\n');
