%% Run Permutation on two vector with equal size
% - Input1: first vector
% - Input2: second vector
% - CTFUN:  string denoting function used to compute difference between
%           groups (e.g., 'nanmean' or 'nanmedian' or 'nanstd').
% - Mask: Mask the lower triangle of the original correlation matrix that
%         you used to create vector for the input
% - Maskidx: Indices for the vectorized Mask
% - nperm:  number of permutations to run (50000 seems fairly stable)
% - PLOTFIGURES: 1 if you want it to plot figures, else use 0
% - q = FDR Threshold value
%
% OUTPUT
%   Results = Structure with a bunch of output from the analysis
%
% Example usage
%
%   % Look at the simple effect of Control>Patients
%   [Results] = Do_Permutation(Control, Patients, 'nanmean', Mask, Maskidx, 50000, 1, 0.05) 
%
%
function [Results] = rb_Permutation(Input1,Input2,CTFUN,Mask,Maskidx,nperm,PLOTFIGURES,q)

%% Probably built in some sort of check to see if vector dimensions match
if size(Input1,2)~=size(Input2,2)
    error('Input1 and Input2 must have the same size in the column dimension');
end

%%
fprintf('Running conversion of r-values to Fishers Z-values\n');
data_A = zscore(Input1);
data_B = zscore(Input2);
AllData = [data_A;data_B];  

% Make function that will compute Central Tendency Measure
ctfun = @(x,y) eval(sprintf('%s(x)',CTFUN)) - eval(sprintf('%s(y)',CTFUN));  % computes difference on CTFUN

% find difference on whatever measure of central tendency you want in real data
CentralTendencyMeasure = zeros(nperm+1,size(data_A,2)); % pre-allocate in memory
CentralTendencyMeasure(1,:) = ctfun(data_A,data_B);  % compute difference between groups

% Run permutations
for iperm = 1:nperm
    fprintf('Running permutation %d\n',iperm);
    ShuffleRows = randperm(size(AllData,2));  % generate random permutation of row indices
    ShuffledData = AllData(ShuffleRows); % grab and sort the data by ShuffleRows
    data_A_perm = ShuffledData(1:size(data_A,1));  % grab a similar sized matrix as data_A
    data_B_perm = ShuffledData(size(data_B,1)+1:end);  % grab a similar sized matrix as data_B
    CentralTendencyMeasure(iperm+1,:) = ctfun(data_A_perm,data_B_perm);  % compute difference
end % for iperm

% Compute p-value for real data compared to the null distribution generated via permutations
for icorr = 1:size(AllData,2)
    perm_pval(icorr) = sum(abs(CentralTendencyMeasure(:,icorr))>=abs(CentralTendencyMeasure(1,icorr))) / (nperm+1);
end

% Get FDR threshold
 % corrected target threshold
[pID pN] = FDR(perm_pval,q);  % pID and pN are different types of FDR thresholds, so the choice is really up to you...

% Find cells that pass FDR threshold
FDRmask = perm_pval<=pID;
if sum(FDRmask) ==0
    fprintf('No cells pass FDR correction\n');
    cdt = [];
    cells_passingFDRcorr = [];
    % Compute Cohen's d
    nA = size(data_A,1); % number of observations in data_A
    nB = size(data_B,1); % number of observations in data_B
    pooled_variances = ((nA-1)*diag(nancov(data_A)) + (nB-1)*diag(nancov(data_B))) / (nA+nB); % pooled variances
    d_values = CentralTendencyMeasure(1,:) ./ sqrt(pooled_variances)';  % univariate effect sizes
elseif sum(FDRmask)>0
    fprintf('%d cells pass FDR correction\n',sum(FDRmask));
    cells_passingFDRcorr = find(FDRmask);
    % x and y are the row column indices for the cells that pass FDR correction
    [x, y] = ind2sub(1,Maskidx(cells_passingFDRcorr));
    % Get the direct comparisons for the matrix indices that pass FDR
    cdt = [x y];
    data_A_GreaterThan_data_B = CentralTendencyMeasure(1,cells_passingFDRcorr)>0;
    data_A_LessThan_data_B = CentralTendencyMeasure(1,cells_passingFDRcorr)<0;   
    % Compute Cohen's d
    nA = size(data_A,1); % number of observations in data_A
    nB = size(data_B,1); % number of observations in data_B
    pooled_variances = ((nA-1)*diag(nancov(data_A)) + (nB-1)*diag(nancov(data_B))) / (nA+nB); % pooled variances
    d_values = CentralTendencyMeasure(1,:) ./ sqrt(pooled_variances)';  % univariate effect sizes
end % if sum(FDRmask)

%% Fill the Results output structure
Results.CTFUN = CTFUN;  % function used as metric of between-group difference
Results.nperm = nperm;  % number of permutations 
Results.FDR_q = q;  % FDR correction level (the number of false positives, e.g., 0.05 means 5% false positives
Results.FDR_pID = pID;  % standard FDR threshold
Results.FDR_pN = pN;  % non-parametric FDR threshold
% Results.Suprathreshold_RowMask shows which cells in a big long row vector are significant.  
% Can use this to overlay onto the Results.ObservedCentralTendencyMeasure 
% to see what the directionality of the result is.  If a positive number it
% means data_A > data_B, while negative number means data_A < data_B.
Results.FDR_Suprathresh_RowMask = FDRmask;  
Results.FDR_SuprathreshCellIndices = cdt;  % x y cell indices of FDR suprathreshold cells
Results.perm_pval = perm_pval;  % permutation pvalues
Results.ObservedCentralTendencyMeasure = CentralTendencyMeasure(1,:);  % real measures of difference between two groups
Results.PermCentralTendencyMeasure = CentralTendencyMeasure(2:end,:);  % measures of difference between two groups from random permutations.
% Results.SuprathresholdCellMatrix_wDirectionality is a matrix the same
% size as the average correlation matrix, filled with zeros anywhere there
% are cells that don't pass FDR correction.  In the cells that do pass FDR
% correction, I've put either 1 or -1 in them.  A 1 means the result was in
% the direction of data_A > data_B, while -1 means the result was in the
% direction of data_A < data_B.  
Results.SuprathresholdCellMatrix_wDirectionality = zeros(size(cdt));
for icell = 1:length(cells_passingFDRcorr)
    if data_A_GreaterThan_data_B(icell) ==1
        Results.SuprathresholdCellMatrix_wDirectionality(cdt(icell,1),cdt(icell,2)) = 1;
    elseif data_A_LessThan_data_B(icell)==1
        Results.SuprathresholdCellMatrix_wDirectionality(cdt(icell,1),cdt(icell,2)) = -1;
    end % if
end % for icell
% Univariate Effect Sizes
Results.Cohens_d = d_values;  % Cohen's d for all cells tested
Results.Cohens_d_FDR_SuprathreshCells = d_values(cells_passingFDRcorr);  % Cohen's d for all cells that pass FDR correction
% Results.SuprathresholdCellMatrix_wEffectSizes is a matrix the same
% size as the average correlation matrix, filled with zeros anywhere there
% are cells that don't pass FDR correction.  In the cells that do pass FDR
% correction, I've put in the effect sizes
Results.SuprathresholdCellMatrix_wEffectSizes = zeros(size(cdt));
for icell = 1:length(cells_passingFDRcorr)
    if data_A_GreaterThan_data_B(icell) ==1
        Results.SuprathresholdCellMatrix_wEffectSizes(cdt(icell,1),cdt(icell,2)) = Results.Cohens_d_FDR_SuprathreshCells(icell);
    elseif data_A_LessThan_data_B(icell)==1
        Results.SuprathresholdCellMatrix_wEffectSizes(cdt(icell,1),cdt(icell,2)) = Results.Cohens_d_FDR_SuprathreshCells(icell);
    end % if
end % for icell

%% Plot some figures to visualize main results
% Plot null distribution from permutations
if PLOTFIGURES
    figure;
    nbins = min(CentralTendencyMeasure):.01:max(CentralTendencyMeasure);
    hist(CentralTendencyMeasure(2:end,:),nbins);
    xlabel(sprintf('%s Difference',CTFUN));
    title('Null Distribution');
    
    % Plot matrix where suprathreshold cells are filled in with their effect size
    figure;
    imagesc(Results.SuprathresholdCellMatrix_wEffectSizes); colorbar;
    title('FDR Suprathreshold cells filled with Cohens d values');
end % if PLOTFIGURES
end % function [Results] = Do_Permutation(Input1, Input2, Mask, Maskidx)
