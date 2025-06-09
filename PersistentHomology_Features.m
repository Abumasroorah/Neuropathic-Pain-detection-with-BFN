clear; 
close all;
clc;

% Define data directory and sampling frequency
dataDir = fullfile(pwd, 'NeuropathicPainEEG');
Fs = 256;            % Sampling rate
L = 3 * Fs;          % 2 seconds = 512 samples

% Initialize output
DATA = cell(49,1);

% Loop through all subjects
for j = 1:size(DATA)
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
save('PainData_3s.mat', 'DATA', '-v7.3');
fprintf('Segmentation complete. Data saved as PainData_3s.mat\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting the BFN fo the 3s_epoched data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear;
load BFN_Pain_3s;

fs = 256;
tw = -1/4 + 1/fs : 1/fs : 2 - 1/fs;
N=size(DATA,2);
% Adjust frequency range logic:
% First 4 subjects: full frequency range (1–100 Hz)
% Remaining subjects: limited to high gamma (45–100 Hz)
Fc = linspace(1,45,90);

for subj = 1:N
    fprintf('Processing Subject %d...\n', subj);
    A1 = double(cell2mat(DATA(subj)));
    A = permute(A1, [3 2 1]);
    B = A;
  
    tic;
    wPLI = wpli(A, B, tw, Fc);
    filename = sprintf('wPLI%d.mat', subj);
    save(filename, 'wPLI', '-v7.3');
    toc;

    clear A1 A B wPLI;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear;
%% STEP 1: Load all PLV files for the eyes closed condition
num_subjects = 49; 
Data1 = cell(num_subjects, 1);

for i = 1:num_subjects
    filename = sprintf('wPLI%d.mat', i);
    data = load(filename);
    varname = sprintf('wPLI_%d', i);
    Data1{i} = data.(varname);
end

%% STEP 2: Band-wise PLV extraction
bandwPLI = cell(num_subjects, 1);
for i = 1:num_subjects
    bandwPLI{i} = bandxtract(Data1{i});  % assumes bandxtract returns a cell array per band
end

%% STEP 3: OMST thresholding per band for formulate BFN
Band_BFN = cell(num_subjects, 1);
for ii = 1:num_subjects
    aa = bandwPLI{ii};  % aa is a cell array where each cell is a band
    for j = 1:numel(aa)  % for each band
        band_data = squeeze(mean(aa{j}, 3));  % average across third dim
        for jj = 1:size(band_data, 3)
            BFN(:, :, jj) = omst_thresholding(band_data(:, :, jj));
        end
        Dat{j} = BFN;  % store thresholded BFN for this band
    end
    Band_BFN{ii} = Dat;
    clear BFN Dat;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computing persistent homology features
neighborhood_radius = 0.5; % You can adjust this parameter based on your data



for ii=1:size(Band_BFN,1)
    a = Band_BFN{ii,1};
    
    for jj=1:size(a,2)
        aa = squeeze(a{1,jj});
        
        for i=1:size(aa,4)
            PLI =squeeze(mean(aa,3));
            distance_mat = 1-squeeze(PLI(:,:,i));
            simplicial_complex = vietoris_rips(distance_mat, neighborhood_radius);
            per_di = compute_persistence_diagram(simplicial_complex);
            sigma = estimate_sigma(per_di);
            PL(i)=mean(persistence_landscape(per_di));
            Bcurv(i) = mean(betti_curve(per_di));
            Ampl(i) = mean(amplitude(per_di));
            Entr(i) = mean(persistent_entropy(per_di));
            NDP(i) = mean(non_diagonal_points(per_di));
            p_kern(i) = mean(pssk_embedding_cell(per_di,sigma));
            pds_kern(i) = mean(persistence_diagram_kernel_cell(per_di,sigma));
            cl(i) = mean(clustering_coef_wu(squeeze(PLI(:,:,i))));
            [pl(i),eff(i,:),ecc(i,:),rad(i,:),~] = charpath(distance_mat);
%             clear Simplificial_complex; clear per_di;
        end
       P(jj) = {PL};
       B(jj) = {Bcurv};
       P_E(jj) = {Entr};
       A(jj) = {Ampl};
       N(jj) = {NDP};
       Psk(jj) = {p_kern};
       pds_k(jj) = {pds_kern};
       Cl(jj) = {cl};
       Pl(jj) = {pl};
       Eff(jj) = {(mean(eff,2))'};
       Ecc(jj) = {(mean(ecc,2))'};
       Rad(jj) = {(mean(rad,2))'};
       
       clear PL; clear Bcurv; clear Ampl; clear Entr; clear NDP;clear p_kern;clear pds_kern rad ecc eff pl cl;
    end
    P_Landscape(ii) = {P};
    Bett_curve(ii) = {B};
    Persist_entr(ii) = {P_E};
    Amplitude(ii) = {A};
    Non_dpoint(ii) = {N};
    P_SK(ii) = {Psk};
    P_DS(ii) = {pds_k};
    Clust(ii) = {Cl};
    ChPL(ii) = {Pl};
    EFF(ii) = {Eff};
    ECC(ii) = {Ecc};
    RAD(ii) = {Rad};
    clear P; clear B; clear P_E; clear A; clear N;clear Pk; clear pds_k Cl Pl Eff Ecc Rad;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Class indices
High_NP_idx     = [1 4 7 10 15 18 19 21 23 24 28 32 33 35];
Moderate_NP_idx = [2 5 9 12 16 17 20 22 29 34 36];
Low_NP_idx      = [3 6 8 11 13 14 25 26 27 30 31];
NC_idx          = 37:49;
Bett_curve_fEC1 = Bett_curve;
Bett_curve_fEC2 = Amplitude;
Bett_curve_fEC3 = P_Landscape;
Bett_curve_fEC4 = Non_dpoint;
Bett_curve_fEC5 = Persist_entr;
Bett_curve_fEC6 = Clust;
Bett_curve_fEC7 = ChPL;
Bett_curve_fEC8 = gEff;
Bett_curve_fEC9 = ECC;
Bett_curve_fEC10 = RAD;



% Initialize
Feature_matrix1 = [];
Feature_matrix2 = [];
Feature_matrix3 = [];
Feature_matrix4 = [];
Feature_matrix5 = [];
Feature_matrix6 = [];
Feature_matrix7 = [];
Feature_matrix8 = [];
Feature_matrix9 = [];
Feature_matrix10 = [];
Label_vector = [];

for subj = 1:49
    % Assign class label
    if ismember(subj, High_NP_idx)
        class_label = 1;
    elseif ismember(subj, Moderate_NP_idx)
        class_label = 2;
    elseif ismember(subj, Low_NP_idx)
        class_label = 3;
    elseif ismember(subj, NC_idx)
        class_label = 4;
    else
        continue;
    end

    % Process all 5 feature sets
    for i = 1:10
        % Get the right Bett_curve
        eval(['feature_cells = Bett_curve_fEC' num2str(i) '{subj};']);

        n = length(feature_cells{1}); % Number of samples for this subject

        % Form [n × 5] matrix
        sample_matrix = zeros(n, 10);
        for f = 1:5
            sample_matrix(:, f) = feature_cells{f}(:);
        end

        % Append to corresponding Feature_matrix
        eval(['Feature_matrix' num2str(i) ' = [Feature_matrix' num2str(i) '; sample_matrix];']);
    end

    % Add labels
    Label_vector = [Label_vector; class_label * ones(n, 1)];
end

% Final concatenation: per-frequency-band matrices
Feat_matrix_delt = [Feature_matrix1(:, 1), Feature_matrix2(:, 1), Feature_matrix3(:, 1), Feature_matrix4(:, 1), Feature_matrix5(:, 1), Feature_matrix6(:, 1), Feature_matrix7(:, 1), Feature_matrix8(:, 1), Feature_matrix9(:, 1), Feature_matrix10(:, 1)];
Feat_matrix_thet = [Feature_matrix1(:, 2), Feature_matrix2(:, 2), Feature_matrix3(:, 2), Feature_matrix4(:, 2), Feature_matrix5(:, 2), Feature_matrix6(:, 2), Feature_matrix7(:, 2), Feature_matrix8(:, 2), Feature_matrix9(:, 2), Feature_matrix10(:, 2)];
Feat_matrix_alp  = [Feature_matrix1(:, 3), Feature_matrix2(:, 3), Feature_matrix3(:, 3), Feature_matrix4(:, 3), Feature_matrix5(:, 3), Feature_matrix6(:, 3), Feature_matrix7(:, 3), Feature_matrix8(:, 3), Feature_matrix9(:, 3), Feature_matrix10(:, 3)];
Feat_matrix_bet  = [Feature_matrix1(:, 4), Feature_matrix2(:, 4), Feature_matrix3(:, 4), Feature_matrix4(:, 4), Feature_matrix5(:, 4), Feature_matrix6(:, 4), Feature_matrix7(:, 4), Feature_matrix8(:, 4), Feature_matrix9(:, 4), Feature_matrix10(:, 4)];
Feat_matrix_gam  = [Feature_matrix1(:, 5), Feature_matrix2(:, 5), Feature_matrix3(:, 5), Feature_matrix4(:, 5), Feature_matrix5(:, 5), Feature_matrix6(:, 5), Feature_matrix7(:, 5), Feature_matrix8(:, 5), Feature_matrix9(:, 5), Feature_matrix10(:, 5)];


filename = 'New_featureEC2';
save(filename, 'Feat_matrix_delt', 'Feat_matrix_thet', 'Feat_matrix_alp', 'Feat_matrix_bet', 'Feat_matrix_gam', 'Label_vector');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all pairs of features extractions

% clear; clc;
% load New_featureEC2
% % Assuming these variables are already loaded
% % Feat_matrix_delt: N x d feature matrix
% % Label_vector:     N x 1 label vector with values 1, 2, 3, or 4
% 
% % Define class labels and names
% class_labels = [1, 2, 3, 4];
% class_names = {'HP', 'MP', 'LP', 'NP'};
% 
% % Generate all binary combinations (2 classes out of 4)
% binary_combos = nchoosek(class_labels, 2);
% 
% % Generate all 3-class combinations (3 classes out of 4)
% three_class_combos = nchoosek(class_labels, 3);
% 
% % Store binary combinations
% fprintf('--- Binary Class Combinations ---\n');
% for i = 1:size(binary_combos, 1)
%     cls = binary_combos(i, :);
%     idx = ismember(Label_vector, cls);
%     X_bin = Feat_matrix_delt(idx, :);
%     y_bin = Label_vector(idx);
% 
%     % Optional: Relabel to 0 and 1 for SVM if needed
%     y_bin_relabel = y_bin == cls(2);  % Assign class 2 as label '1', class 1 as '0'
% 
%     % Save to struct or workspace
%     varname = sprintf('X_%s_vs_%s', class_names{cls(1)}, class_names{cls(2)});
%     assignin('base', varname, X_bin);
%     varname_y = sprintf('y_%s_vs_%s', class_names{cls(1)}, class_names{cls(2)});
%     assignin('base', varname_y, y_bin_relabel);
% 
%     fprintf('Created: %s, %s\n', varname, varname_y);
% end
% 
% % Store 3-class combinations
% fprintf('\n--- Three-Class Combinations ---\n');
% for i = 1:size(three_class_combos, 1)
%     cls = three_class_combos(i, :);
%     idx = ismember(Label_vector, cls);
%     X_multi = Feat_matrix_delt(idx, :);
%     y_multi = Label_vector(idx);  % Keep original labels (1/2/3/4) for multiclass
% 
%     % Save to struct or workspace
%     varname = sprintf('X_%s_%s_%s', class_names{cls(1)}, class_names{cls(2)}, class_names{cls(3)});
%     assignin('base', varname, X_multi);
%     varname_y = sprintf('y_%s_%s_%s', class_names{cls(1)}, class_names{cls(2)}, class_names{cls(3)});
%     assignin('base', varname_y, y_multi);
% 
%     fprintf('Created: %s, %s\n', varname, varname_y);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the extracted Features and labels as 
clear; clc;
load New_featureEC2

% Assuming these variables are already loaded:
% Feat_matrix_delt: N x d feature matrix
% Label_vector:     N x 1 label vector with values 1, 2, 3, or 4

% Define class labels and names
class_labels = [1, 2, 3, 4];
class_names = {'HP', 'MP', 'LP', 'NP'};

% Generate all binary combinations (2 classes out of 4)
binary_combos = nchoosek(class_labels, 2);

% Generate all 3-class combinations (3 classes out of 4)
three_class_combos = nchoosek(class_labels, 3);

% Initialize feature structure
featureEC = struct();

% Store binary combinations
fprintf('--- Binary Class Combinations ---\n');
for i = 1:size(binary_combos, 1)
    cls = binary_combos(i, :);
    idx = ismember(Label_vector, cls);
    X_bin = Feat_matrix_delt(idx, :);
    y_bin = Label_vector(idx);
    y_bin_relabel = y_bin == cls(2);  % Assign class 2 as label '1', class 1 as '0'

    % Construct field names
    field_x = sprintf('X_%s_vs_%s', class_names{cls(1)}, class_names{cls(2)});
    field_y = sprintf('y_%s_vs_%s', class_names{cls(1)}, class_names{cls(2)});

    % Store in struct
    featureEC.(field_x) = X_bin;
    featureEC.(field_y) = y_bin_relabel;

    fprintf('Created: %s, %s\n', field_x, field_y);
end

% Store 3-class combinations
fprintf('\n--- Three-Class Combinations ---\n');
for i = 1:size(three_class_combos, 1)
    cls = three_class_combos(i, :);
    idx = ismember(Label_vector, cls);
    X_multi = Feat_matrix_delt(idx, :);
    y_multi = Label_vector(idx);  % Keep original labels

    % Construct field names
    field_x = sprintf('X_%s_%s_%s', class_names{cls(1)}, class_names{cls(2)}, class_names{cls(3)});
    field_y = sprintf('y_%s_%s_%s', class_names{cls(1)}, class_names{cls(2)}, class_names{cls(3)});

    % Store in struct
    featureEC.(field_x) = X_multi;
    featureEC.(field_y) = y_multi;

    fprintf('Created: %s, %s\n', field_x, field_y);
end

% Save the full structure
save('featureEC.mat', 'featureEC');
fprintf('\nSaved all features to featureEC.mat\n');









































