function BFN_hec = c2array(BFN_HEC)

% Converting cell of cell to 4D (L by n by n by N) where L is the bands, n
% by n is the adjacency matrix and N is the total number of samples across 
% all subjects.

numSubjects = numel(BFN_HEC);
numLayers = length(BFN_HEC{1,1});
a= BFN_HEC{1,1};
aa=a{1,1};

matrixSize = size(aa,1);

% Determine total number of samples (N)
N = 0;
for subj = 1:numSubjects
    for layer = 1:numLayers
        N = N + size(BFN_HEC{subj}{layer}, 3);
    end
end

% Initialize output matrix
BFN_hec = zeros(numLayers, matrixSize, matrixSize, N);

% Fill BFN_hec
idx = 1;
for subj = 1:numSubjects
    for layer = 1:numLayers
        tempMat = BFN_HEC{subj}{layer};  % size: 22 × 22 × n
        n = size(tempMat, 3);
        BFN_hec(layer, :, :, idx:idx+n-1) = tempMat;
    end
    idx = idx + n;  % Move index forward by n
end




