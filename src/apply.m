
%% *****  LOAD DATA  ******************************************************

fprintf(1,'\n*****  LOAD APPLICATION DATA \n\n')

% set data file to read
dft  = 'MODEL';
mdl = input(['->  What dataset would you like to apply the model to? \n', ...
              '    Give full name of datafile (dft = MODEL): \n'],'s');
if isempty(mdl); mdl = dft; end

% read file to load dataset
if strcmp(mdl(end-3:end),'.csv') || ...  % load data from text file
   strcmp(mdl(end-3:end),'.txt') || ...
   strcmp(mdl(end-3:end),'.dat')
    DATA   = readtable(mdl);
    PRJCT  = DATA.Properties.VariableNames(1    ); PRJCT = PRJCT{:};
    VNAMES = DATA.Properties.VariableNames(2:end);
    SNAMES = DATA{:    ,1};
    Xa      = DATA{:,2:end};
    if any(sum(Xa,2) > 2)  % detect that chemical data is in wt %, convert to wt fraction
        Xa = Xa./100;
        if exist('Ft','var'); Ft = Ft./100; end
        if exist('Xt','var'); Ft = Ft./100; end
    end
    vstype = 'sct';  % scatter plots
elseif strcmp(mdl(end-3:end),'.tif') || ...  % load image data from tiff file
       strcmp(mdl(end-3:end),'tiff')
    T   = Tiff(mdl,'r');
    IMG = read(T);
    if any(size(IMG(:)) > 1e5)  % coarsen image to avoid working with oversized data
        dft = 1;
        crs = input(['\nImage has dimensions of ',int2str(size(IMG,1)),' x ',int2str(size(IMG,2)),' x ' ,int2str(size(IMG,3)),'; apply coarsening factor (coarsen>1, dft=1):\n']);
        if isempty(crs); crs = dft; end
        IMG = IMG(1:crs:end,1:crs:end,:);
        for i = 1:length(size(IMG,3))
           IMG(:,:,i) = imnlmfilt(IMG(:,:,i));
        end
    end
    [mx,my,n] = size(IMG);
    m   = mx*my;
    Xa  = reshape(IMG,m,size(IMG,3));
    Xa  = double(Xa);
    VNAMES = cell(1,n); for i=1:n; VNAMES(i) = {['BND ',int2str(i)]}; end
    SNAMES = {};
    vstype = 'img';  % image plots
else % assume data is provided as Matlab file
    load(mdl);
    if any(sum(Xa,2) > 2)  % detect that chemical data is in wt %, convert to wt fraction
        Xa = X./100;
        if exist('Ft','var'); Ft = Ft./100; end
        if exist('Xt','var'); Xt = Xt./100; end
    end
    vstype = 'sct';  % scatter plots
end
if ~exist('PRJCT','var'); filesplt = split(mdl,'.'); PRJCT = filesplt{1}; end


%% *****  APPLY MODEL  ****************************************************

DGNa = DGN;
[DGNa.m,DGNa.n]  = size(X);
DGN.Ii = (1:DGNa.m).';
DGN.Ir = [];
DGN.fn = 1;

[DGNa.mx,DGNa.my,DGNa.n] = size(IMG);
DGNa.m  = DGNa.mx*DGNa.my;
DGNa.fn = DGNa.fn;

DGNa.fh(DGNa.fn) = visualise(DGNa.fn,{Xa},{'data'},'Application Data',DGNa,VNAMES,vstype); DGNa.fn = DGNa.fn+1;
if strcmp(vstype,'img')
    DGNa.fh(DGNa.fn) = visualise(DGNa.fn,{Xa},{'data'},'Application Data',DGNa,VNAMES,'rgb'); DGNa.fn = DGNa.fn+1;
end

% scale application data to relative intensities/concentrations
Xa = Xa./sum(Xa,2);

% predict PC coordinates by least-squares
Aap = (Xa - DGN.meanX)/DGN.PC(1:DGN.p-1,:);
Xap = Aap*DGN.PC(1:DGN.p-1,:) + DGN.meanX;

% predict EM abundances by least-squares
if iem;  Aai = Xa/Fi;  end
if eem;  Aae = Xa/Fe;  end

% predict cluster membership by nearest distance
[~,Ic] = pdist2(Fc,Xa,'euclidean','Smallest',1);  Ic = Ic.';

if strcmp(vstype,'img')
            DGNa.fh(DGN.fn)  = visualise(DGNa.fn,{Xap},{},'Reduced Data',DGNa,VNAMES,'rgb'); DGNa.fn = DGNa.fn+1;
    if iem; DGNa.fh(DGNa.fn) = visualise(DGNa.fn,{Aai},{},['Map of ',num2str(DGNa.p),' internal EM proportions;'],DGNa,VNAMES,'rgb');  DGNa.fn = DGNa.fn+1; end
    if eem; DGNa.fh(DGNa.fn) = visualise(DGNa.fn,{Aae},{},['Map of ',num2str(DGNa.p),' external EM proportions;'],DGNa,VNAMES,'rgb');  DGNa.fn = DGNa.fn+1; end
    if clu; DGNa.fh(DGNa.fn) = visualise(DGNa.fn,{Ic} ,{},['Map of ',num2str(DGNa.c),' clusters;'],DGNa,VNAMES,'img');  DGNa.fn = DGNa.fn+1; end
end
