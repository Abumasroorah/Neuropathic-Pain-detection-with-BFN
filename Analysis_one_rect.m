clear;
clc;

%% STEP 1: Load all wPLI files
num_subjects = 49; 
% for all the NP subjects and the 13 normal control 
% subjects (eyes closed condition). Similar processing is required for eyes
% opened codition, otherwise, both condition could be processed in parallel.


Data1 = cell(num_subjects, 1);

for i = 1:num_subjects
    filename = sprintf('wPLI_%d.mat', i);
    data = load(filename);
    varname = sprintf('wPLI_%d', i);
    Data1{i} = data.(varname);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot of Time frequency maps for a channel pair and a subject per group
%%%% Drawing the time frequency wPLI of some subjects;
aa1= squeeze(wPLI_1(3,7,:,1:44));
aa2= squeeze(wPLI_2(3,7,:,1:33));
aa3= squeeze(wPLI_3(3,7,:,1:50));
aa4= squeeze(wPLI_38(3,7,:,:));

figure
subplot(2,2,1)
imagesc(abs(a1))
set(gca, 'YDir', 'normal');
subplot(2,2,2)
imagesc(abs(a2))
set(gca, 'YDir', 'normal');
subplot(2,2,3)
imagesc(abs(a3))
set(gca, 'YDir', 'normal');
subplot(2,2,4)
imagesc(abs(a4))
set(gca, 'YDir', 'normal');


%% High_NP, Moderate_NP, Low_NP and NC indices (1-based), the indices could
% be obtained automatically from the file included in the data repository.
High_NP_idx     = [1 4 7 10 15 18 19 21 23 24 28 32 33 35];
Moderate_NP_idx = [2 5 9 12 16 17 20 22 29 34 36];
Low_NP_idx      = [3 6 8 11 13 14 25 26 27 30 31];
NC_idx = 37:49; 

HighNP_wPLI = Data(High_NP_idx);
ModNP_wPLI = Data(Moderate_NP_idx);
LowNP_wPLI = Data(Low_NP_idx);
NC_wPLI = Data(NC_idx);


%% STEP 2: Band-wise PLV extraction
bandwPLI = cell(num_subjects, 1);
for i = 1:num_subjects
    bandwPLI{i} = bandxtract(Data1{i});  % assumes bandxtract returns a cell array per band
end

%% STEP 3: OMST thresholding per band
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

%% STEP 4: Group assignment
% High_NP, Moderate_NP, Low_NP and NC indices (1-based), the indices could
% be obtained automatically from the file included in the data repository.
% High_NP_idx     = [1 4 7 10 15 18 19 21 23 24 28 32 33 35];
% Moderate_NP_idx = [2 5 9 12 16 17 20 22 29 34 36];
% Low_NP_idx      = [3 6 8 11 13 14 25 26 27 30 31];
% NC_idx = 37:49; 

BFN_HEC = Band_BFN(High_NP_idx);
BFN_MEC = Band_BFN(Moderate_NP_idx);
BFN_LEC = Band_BFN(Low_NP_idx);
BFN_NEC = Band_BFN(NC_idx);

%% STEP 5: Convert to standard array format for downstream analysis/visualization
BFN_hec = c2array(BFN_HEC);
BFN_mec = c2array(BFN_MEC);
BFN_lec = c2array(BFN_LEC);
BFN_nec = c2array(BFN_NEC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  For average BFN plots of all the groups.
BandwPLIEC = [];
for ii=1:size(bandwPLI,2)
    a = bandwPLI{1,ii};
    wPLI=[];
    for jj=1:size(a,2)
        aa = a{1,jj};
        PLI =squeeze(mean(aa,3));
        wPLI(jj,:,:,:)=PLI;
       
    end
    clear aa; clear a; clear PLI;
    BandwPLIEC(ii) = {wPLI};
    clear wPLI;
end
    



BandBFN_HEC = BandwPLIEC(1,High_NP_idx); %High_NP group
Band_BFN_HEC = [];
% Loop and concatenate along the 4th dimension
for i = 1:length(BandBFN_HEC)
    Band_BFN_HEC = cat(4, Band_BFN_HEC, BandBFN_HEC{i});
end


BandBFN_MEC = BandwPLIEC(1,Moderate_NP_idx); %Moderate_NP group
Band_BFN_MEC = [];
% Loop and concatenate along the 4th dimension
for i = 1:length(BandBFN_MEC)
    Band_BFN_MEC = cat(4, Band_BFN_MEC, BandBFN_MEC{i});
end


BandBFN_LEC = BandwPLIEC(1, Low_NP_idx); %Low_NP group
Band_BFN_LEC = [];
% Loop and concatenate along the 4th dimension
for i = 1:length(BandBFN_LEC)
    Band_BFN_LEC = cat(4, Band_BFN_LEC, BandBFN_LEC{i});
end



Av_HEC = squeeze(mean(Band_BFN_HEC,4));
Av_HEC = Av_HEC(:, 1:18, 1:19); 
Av_MEC = squeeze(mean(Band_BFN_MEC,4));
Av_MEC = Av_MEC(:,1:18, 1:19); 

Av_LEC = squeeze(mean(Band_BFN_LEC,4));
Av_LEC = Av_LEC(:,1:18, 1:19); 





for ii=1:size(Av_HEC,1)
    HEC(ii,:,:) = omst_thresholding(squeeze(Av_HEC(ii,:,:)));
	MEC(ii,:,:) = omst_thresholding(squeeze(Av_MEC(ii,:,:)));
    LEC(ii,:,:) = omst_thresholding(squeeze(Av_LEC(ii,:,:)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%Here is the plot for the High NP group at Eyes Closed Condition, 
% Plot for other conditions could be generated in a similar manner.

for kk=1:size(HEC,1)
    for ii=1:size(HEC,2)
        for jj=1:size(HEC,3)
            if(HEC(kk,ii,jj)>0)
                HEC_MRPs(kk,ii,jj) = 1;
            else
                HEC_MRPs(kk,ii,jj) = 0;
            end
        end
    end
end

for ii=1:size(HEC_MRPs,1)
    figure(ii)
    config.layout='Pain19channn.lay';
    lay=Layout(config);
    plot_topo2(lay,(squeeze(HEC_MRPs(ii,:,:))),1);
    colormap winter
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculating the groupwise iterative meanimum dominating sets

%% Loading Data1
clear;
close all
clc;
load Data1

for ii=1:size(Data1,1)
    a = Data1{ii, 1};
    
%     for jj=1:size(a,2)
%         aa = a{1,jj};
        wPLI =squeeze(mean(a,3));
       
%     end
    clear aa; clear a; clear PLI;
    BandwPLIEC(ii) = {wPLI};
    clear wPLI;
end
   



for ii=1:size(BandwPLIEC,2)
    a = BandwPLIEC{1, ii};
    
%     for jj=1:size(a,2)
%         aa = a{1,jj};
        PLI =squeeze(mean(a,3));
%         PLI = squeeze(mean(PLI,3));
        wPLI=PLI;
       
%     end
    clear a; clear PLI;
    Band_wPLIEC(ii,:,:) = wPLI;

    clear wPLI;
end


High_NP_idx     = [1 4 7 10 15 18 19 21 23 24 28 32 33 35];
Moderate_NP_idx = [2 5 9 12 16 17 20 22 29 34 36];
Low_NP_idx      = [3 6 8 11 13 14 25 26 27 30 31];
NC_idx = 37:49; 


BandBFNEC_HighEC = Band_wPLIEC(High_NP_idx,:,:);
BandBFNEC_ModEC = Band_wPLIEC(High_NP_idx,:,:);
BandBFNEC_LowEC = Band_wPLIEC(High_NP_idx,:,:);



%% Separating into groups
Av_HEC = squeeze(BandBFNEC_HighEC(:,1:18, 1:19));
Av_MEC = squeeze(BandBFNEC_ModEC(:,1:18, 1:19));
Av_LEC = squeeze(BandBFNEC_LowEC(:,1:18, 1:19));

% Av_NEC = squeeze(mean(BandBFNEC_NCEC,3));



% Av_HEO = squeeze(BandBFNEC_HighEO(:,1:18, 1:19));
% Av_MEO = squeeze(BandBFNEC_ModEO(:,1:18, 1:19));
% Av_LEO = squeeze(BandBFNEC_LowEO(:,1:18, 1:19));

% Av_NEO = squeeze(mean(BandBFNEC_NCEO,3));

for ii=1:size(Av_HEC,1)
    HEC(ii,:,:) = omst_thresholding(squeeze(Av_HEC(ii,:,:)));
end

for ii=1:size(Av_MEC,1)
    MEC(ii,:,:) = omst_thresholding(squeeze(Av_MEC(ii,:,:)));
end

for ii=1:size(Av_LEC,1)
    LEC(ii,:,:) = omst_thresholding(squeeze(Av_LEC(ii,:,:)));
end




for ii=1:size(HEC,1)
    MDS_HighEC = iterativeMinDomSet(squeeze(HEC(ii,:,:)), 5);
    MDS_HEC(ii) = {fiveconf(MDS_HighEC)};
end


MDS_HEC = fiveconf(MDS_HEC);




for ii=1:size(HEC,1)
    a=(squeeze(HEC(ii,:,:)));
    B_grHC{ii} = a;
end


HEC_conncomp = find_frequent_connected_subgraph(B_grHC,.5);

Freq_HEC=sum (HEC_conncomp);

dexHEC = find(Freq_HEC ~= 0);
MDS_sig_HEC = intersect(unique(MDS_HEC), find(Freq_HEC ~= 0));


% highlightedIndices = [1, 2, 3, 4]; % Example indices
MRPs = zeros(19,19);
% MDS_sig_HEC = [5, 10];
    
% A = MDS_sig_HEC;
% highlightedIndices = (A{1,1}); % Example indices
    
figure
config.layout='Pain19channn.lay';
lay=Layout(config);
plot_topo_abdul(lay,1,MDS_sig_HEC, 'm', 120);
colormap winter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:size(MEC,1)
    MDS_ModEC = iterativeMinDomSet(squeeze(MEC(ii,:,:)), 5);
    MDS_MEC(ii) = {fiveconf(MDS_ModEC)};
end


MDS_MEC = fiveconf(MDS_MEC);




for ii=1:size(MEC,1)
    a=(squeeze(MEC(ii,:,:)));
    B_grMC{ii} = a;
end


MEC_conncomp = find_frequent_connected_subgraph(B_grMC,.5);

Freq_MEC=sum (MEC_conncomp);

dexMEC = find(Freq_MEC ~= 0);
MDS_sig_MEC = intersect(unique(MDS_MEC), find(Freq_MEC ~= 0));


% highlightedIndices = [1, 2, 3, 4]; % Example indices
MRPs = zeros(19,19);
% MDS_sig_HEC = [5, 10];
    
% A = MDS_sig_HEC;
% highlightedIndices = (A{1,1}); % Example indices
    
figure
config.layout='Pain19channn.lay';
lay=Layout(config);
plot_topo_abdul(lay,1,MDS_sig_MEC, 'm', 120);
colormap winter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for ii=1:size(LEC,1)
    MDS_LowEC = iterativeMinDomSet(squeeze(LEC(ii,:,:)), 10);
    MDS_LEC(ii) = {fiveconf(MDS_LowEC)};
end


MDS_LEC = fiveconf(MDS_LEC);




for ii=1:size(LEC,1)
    a=(squeeze(LEC(ii,:,:)));
    B_grLC{ii} = a;
end


LEC_conncomp = find_frequent_connected_subgraph(B_grLC,.5);

Freq_LEC=sum (LEC_conncomp);

dexLEC = find(Freq_LEC ~= 0);
MDS_sig_LEC = intersect(unique(MDS_LEC), find(Freq_LEC ~= 0));


% highlightedIndices = [1, 2, 3, 4]; % Example indices
MRPs = zeros(19,19);
% MDS_sig_HEC = [5, 10];
    
% A = MDS_sig_HEC;
% highlightedIndices = (A{1,1}); % Example indices
    
figure
config.layout='Pain19channn.lay';
lay=Layout(config);
plot_topo_abdul(lay,1,MDS_sig_LEC, 'm', 120);
colormap winter


