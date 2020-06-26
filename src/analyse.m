% analyse:  unmix subroutine to analyse data for principal component decomposition
%
%   [X,DGN] = analyse(X,DGN)
%
%   The routine performs principal component analysis (PCA) to prepare reduction 
%   dimensionality. The PCA employs singular value decomposition on data
%   that is centred (reduced by mean) and weighted by variable-wise
%   variances. Key outputs of PCA, including principal component vectors,
%   cumulative variances, and goodness of fit coefficients are calculated
%   for any applicable (p-1)-dimensional reduced data space and passed back
%   in the diagnostics structure.
%
%   X      : input unprocessed data in FCM space
%   DGN    : input structure containing data diagnostics
%
%   X      : output sum-normalised data in FCM space
%   DGN    : output structure with added PCA and projections diagnostics
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function  [X,DGN] = analyse(X,DGN)

X  =  X./sum(X,2);  % normalise data to unit sum

% store variable-wise descriptive statistics for sum-normalised data
DGN.minX  = min(X);
DGN.maxX  = max(X);
DGN.meanX = mean(X);
DGN.mednX = median(X);
DGN.stdvX = std(X);

% test projection into principal component space for p = 2:n components
[PC_C, PC_A, PC_V] = pca(X,'Algorithm','svd','Centered','on','VariableWeights','variance');
DGN.PC = PC_C.';
DGN.PA = PC_A;
DGN.PV = PC_V;
DGN.CV = cumsum(PC_V(1:end-1)./sum(PC_V(1:end-1)));
fp     = zeros(DGN.n-1,1);
for p  = 2:DGN.n
    Fp = DGN.PC(1:p-1,:);
    Ap = DGN.PA(:,1:p-1);
    Xp = Ap*Fp + mean(X);
    DGN.CD(:,p-1) = 1 - std(Xp-X).^2./std(X).^2;  % get correlation coefficients of data fit
    [~,~,sumD]    = kmeans(X,p);
    fp(p-1)       = mean(sumD).^2;
end  % end for
DGN.CC = cumsum(fp./sum(fp));

DGN.CS = DGN.CV.*mean(DGN.CD).'.*DGN.CC;


end  % end function