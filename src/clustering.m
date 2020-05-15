% cluster:  unmix subroutine to analyse clustering of data
%
%   [Ic,Fc,Fcp,DGN] = clustering(A,DGN)
%
%   The routine performs clustering analysis (PCA) in reduced principal
%   component space (RCP). The user can select from different clustering
%   methods, including k-mean, fuzzy c-means, and hierarchical clustering
%   with various distance measures (see help to 'linkage' for options). The
%   routine returns cluster membership indices for each sample and cluster
%   centroid compositions in RCP and full measurement space.
%
%   Ap     : input abundances in RCP space
%   DGN    : input structure containing data diagnostics
%
%   Ic     : output cluster membership indices for samples
%   Fc     : output cluster centroid compositions in FMC space
%   Fcp    : output cluster centroid compositions in RCP space
%   DGN    : output diagnostics structure with added clustering method
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function    [Ic,Fc,Fcp,DGN] = clustering(Ap,DGN)

% select method for clustering analysis
dft = 'k-means';
mth = input(['->  Select method for clustering analysis: \n' ...
             '    k-means clustering: enter k, kmeans, k-means (dft) \n' ...
             '    fuzzy c-means clustering: enter f, fuzzy, fcm \n' ...
             '    hierarchical clustering including linkage method: enter h-[linkage], hierarchical-[linkage], \n' ...
             '    for example, h-single, h-average, h-ward, h-centroid \n'],'s');
if isempty(mth); mth = dft; end

Fcp = zeros(DGN.p,DGN.p-1);

switch mth(1)
    
    case 'f'
        
        mth = 'fuzzy-c-means';

        % perform fuzzy clustering analysis (requires fuzzy toolbox)
        [Fcp,Mc]   = fcm(Ap(DGN.Ii,:),DGN.p,[2,1e3,1e-6,0]);  % clustering in RPC
        Ic = zeros(size(DGN.Ii));
        for c = 1:DGN.p
            Ic(Mc(c,:).'>0.5) = c;
        end

    case 'k'
        
        mth = 'k-means';
                
        % perform k-means clustering analysis 
        [Ic,Fcp] = kmeans(Ap(DGN.Ii,:),DGN.p);

    case 'h'
        
        % performs hierarchical clustering analysis
        link    = split(mth,'-');
        link(1) = {'linkage'};
        if length(link)<2; link = {'linkage','ward'}; end
        mth = ['hierarchical-',link(2)];
        Ic = clusterdata(Ap(DGN.Ii,:),link{:},'maxclust',DGN.p);
        for c = 1:DGN.p
            Fcp(c,:) = mean(Ap(DGN.Ii(Ic==c),:));
        end
        
end

Fc = Fcp*DGN.PC(1:DGN.p-1,:) + DGN.meanX;  % transform back to FMC

DGN.cluster_method = mth;

end  % end function