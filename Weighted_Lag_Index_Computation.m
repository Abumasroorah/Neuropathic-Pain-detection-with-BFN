% Here is the main code for the computation of the wPLI matrices for all
% the subjects in eyes closed condition. Same process should be done for
% eyes open condition.

clear; 
close all; 
clc;

load PainData_6s;
fs = 256;
tw = -1/4 + 1/fs : 1/fs : 2 - 1/fs;
N=size(DATA2,2);
% Adjust frequency range logic:
% First 4 subjects: full frequency range (1–100 Hz)
% Remaining subjects: limited to high gamma (45–100 Hz)
Fc = linspace(1,45,90);

for subj = 1:N
    fprintf('Processing Subject %d...\n', subj);
    A1 = double(cell2mat(DATA2(subj)));
    A = permute(A1, [3 2 1]);
    B = A;
    
    
    
  
    tic;
    wPLI = wpli(A, B, tw, Fc);
    filename = sprintf('wPLI_%d.mat', subj);
%     save(filename, 'G:\Pain', 'PLV', '-v7.3');
    filepath = fullfile('path', filename);
    save(filepath, 'wPLI', '-v7.3');
    toc;

    clear A1 A B wPLI;
end
