% reduce:  unmix subroutine to reduce dimensionality of data set
%
%   [Xp,Ap,DGN] = reduce(X,DGN,VNAMES)
%
%   The routine uses the results of principal component analysis to reduce
%   the input data to a (p-1)-dimensional principal component space.
%   User-specified tolerances are used to decide the optimal number of
%   dimensions. Optimality is measured as a trade-off between number of
%   components (fewer components preferred) and the data variance retained
%   in the reduced space (higher explained variance preferred), and the
%   goodness of fit between data projected to reduced principal component
%   (RCP) compared to data in full measurement space (FCM) (better fit is 
%   preferred).
%
%   X      : input sum-normalised data in FCM space
%   DGN    : input structure containing data diagnostics
%   VNAMES : input cell array with stored variable names for plot labels
%
%   Xp     : output data projected to (p-1)-dimensional RCP space
%   Ap     : output PC abundances for projected data in RCP space
%   DGN    : output structure with added PCA and projections diagnostics
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function  [Xp,Fp,Ap,DGN] = reduce(X,DGN,VNAMES)

% set tolerances for EM selections
dft = [0.9,2.0];
tol = input(['->  Adjust endmember selection tolerances as list [CStol,ORtol] \n' ...
             '    CStol: tolerance for selection coefficients (dft = ',num2str(dft(1)),') \n' ...
             '    ORtol: tolerance for outlier data removal   (dft = ',num2str(dft(2)),') \n' ]);
if isempty(tol); tol = dft; end
CStol  = tol(1);
ORtol  = tol(2);

% select number of EMs 'p'
p      = find(    DGN.CS > CStol,1) + 1;
DGN.p  = p;

% truncate to p-1 principal components and initial outlier removal
repeat = 1;
while repeat
    % truncate to p-1 principal components
    Fp = DGN.PC(1:DGN.p-1,:);
    Ap = DGN.PA(:,1:DGN.p-1);
    Xp = Ap*Fp + mean(X);

    % remove outlier points with bad fit to p-EM model
    ir         = find(any(abs(Xp-X)./std(Xp) > ORtol,2));
    DGN.Ir     = ir;
    DGN.Ii     = (1:DGN.m).';
    DGN.Ii(ir) = [];
    DGN.rm     = length(DGN.Ir);
    
    % visualise principle components and fitted data for p-EM model
    PCNAMES = cell(1,p); for i=1:p; PCNAMES(i) = {['PC',int2str(i)]}; end
    DGN.fh(DGN.fn) = visualise(DGN.fn,{Ap},{'reduced data'},['Reduced data in ',num2str(p),' principal component space;'],DGN,PCNAMES); DGN.fn = DGN.fn+1; 
    
    % plot residuals of data projection
    DGN.fh(DGN.fn) = visualise(DGN.fn,{(Xp-X)./std(Xp)},{'proj. residuals'},['Data projection residuals for ',num2str(DGN.p),' EMs'],DGN,VNAMES); DGN.fn = DGN.fn+1; 

    % plot dimensionality selection metrics
    DGN.fh(DGN.fn) = figure(DGN.fn); clf;  DGN.fn = DGN.fn+1; 
    FS = {'FontSize',14}; MS = {'MarkerSize',8}; LW = {'LineWidth',1.5};
    sgtitle(['Selected ',int2str(p),' endmembers for mixing model'],FS{:})
    subplot(1,3,1)
    plot(2:DGN.n,DGN.CV,'k-',2:DGN.n,DGN.CV,'ko',MS{:},LW{:}); hold on; box on; axis tight;
    plot(p,DGN.CV(p-1),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('data variance',FS{:})
    subplot(1,3,2)
    plot(2:DGN.n,mean(DGN.CD),'k-',2:DGN.n,mean(DGN.CD),'ko',MS{:},LW{:}); hold on; box on; axis tight;
    plot(p,mean(DGN.CD(:,p-1)),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('data fit',FS{:})
    subplot(1,3,3)
    plot([2,DGN.n],[CStol,CStol],'b-',MS{:},LW{:}); hold on; box on; axis tight;
    plot(2:DGN.n,DGN.CS,'k-',2:DGN.n,DGN.CS,'ko',MS{:},LW{:});
    plot(p,DGN.CS(p-1),'ro',MS{:},LW{:});
    set(gca,LW{:});
    title('selection',FS{:})
    
    % decide to procede with selected model or change selection
    dft  = 1;
    prcd = input(['\n->  Proceed with ',num2str(p),'-EM model (1, dft) or adjust number of EM (>1)? \n']);
    if isempty(prcd); prcd = dft; end
    
    if  prcd>1;  p = prcd;   end
    if  prcd==1; repeat = 0; end
    DGN.p  = p;
end  % end while

end  % end function