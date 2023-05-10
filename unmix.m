% unmix:    A tool to interrogate geochemical whole rock analyses through 
%           principal component analysis, clustering analysis and end-member 
%           mixing model fitting.
%
% This script performs principal component analysis, fuzzy cluster extraction, 
% mutually extreme sample identification, and external mixing end-member
% estimation. The routine takes as input a set of geochemical whole rock
% analyses that may contain both major and trace element compositions.
%
% In the first step, the script performs a principal component analysis to 
% selects an optimal number of components to reduce dimensionality and 
% identify outlier samples. The routine projects the data into reduced
% principal component space (RPC) and quantifies goodness of fit coefficients 
% for each chemical variable.
%
% In the second step, the script performs a fuzzy clustering analysis to
% determine how the data clusters in reduced principal component space. The
% routine returns fuzzy cluster centre coordinates in RPC and full
% measurement component (FMC) space, along with cluster membership matrix
% for all samples.
%
% In the third step, the script extracts mutually extreme samples from the
% data set by applying the Vertex Component Analysis (VCA) algorithm developed 
% by [REF] while maximising the volume of the simplex spanned by the
% identified extreme samples in reduced principal component space. The
% routine returns extreme sample indices and coordinates in RPC and FMC
% spaces.
%
% In the fourth step, the script calculates external mixing end-members and
% corresponding abundances that fit the data as projected into RPC space by
% minimising the volume of the data-inclusive simplex spanned by the end-
% members while maintaining non-negative abundances and end-member 
% compositions to a set tolerance. The routine is based on the MINVEST 
% algorithm developed by [REF] and returns end-member coordinates and 
% abundances in RPC and MFC space.
% 
% The script takes as input tabulated data placed in the sub-directory 'data'. 
% Input data should be formatted as text files (*.csv, *.txt, *.dat), with 
% analyses given in [wt %] with one row per sample, one column per chemical 
% variable (major oxide, trace element, etc.), with project tag or dataset 
% name placed in cell [1,1] (first entry on top left), variable names in 
% first row, and sample names/numbers in first column, for example:
%
% EXAMPLE DATA TABLE:
% 
% Example_Data,  VAR_1,   VAR_2,   VAR_3,   ...
% SMP_01,         1.00,    2.00,    3.00,   ...
% SMP_02,         5.00,    4.00,    3.00,   ...
% SMP_03,        10.00,   20.00,   30.00,   ...
% ...              ...      ...      ...    ...
%
% Alternatively, the script takes input files in Matlab format (.mat) with
% data placed in array 'X' of dimensions (# samples x # variables), project 
% tag or dataset name in string 'PRJCT', and sample and variable names in 
% cell arrays SNAMES, and VNAMES, respectively. If using synthetic data 
% generated for benchmarking purposes, use array 'Xt' to provide true data 
% values, 'At' for true mixing proportions, and 'Ft' for true mixing end-
% member compositions.
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


%% *****  PREP WORKSPACE  *************************************************
close all;             % close figures
clear DGN X Xp Ap Fe Fi Fc
warning('off','all');  % turn off warning
addpath ./src ./data   % add paths to source code and data directories

fprintf(1,'\n\n*****  UNMIX - TOOL FOR GEOCHEMICAL AND MULTISPECTRAL DATA ANALYSIS \n\n')


%% *****  LOAD DATA  ******************************************************

fprintf(1,'\n*****  LOAD TRAINING DATA \n\n')

% set data file to read
dft  = 'DATA';
file = input(['->  What dataset would you like to process? \n', ...
              '    Give full name of datafile (dft = DATA): \n'],'s');
if isempty(file); file = dft; end

% read file to load dataset
if exist(file,'var')
    if isfield(DATA,'Properties')
        PRJCT  = DATA.Properties.VariableNames(1); PRJCT = PRJCT{:};
        VNAMES = DATA.Properties.VariableNames(2:end);
        SNAMES = DATA{:,1};
        X      = DATA{:,2:end};
    else
        PRJCT  = DATA.PRJCT;
        VNAMES = DATA.VNAMES;
        SNAMES = DATA.SNAMES;
        X      = DATA.X;
    end
    if any(sum(X,2) > 2)  % detect that chemical data is in wt %, convert to wt fraction
        X = X./100;
        if exist('Ft','var'); Ft = Ft./100; end
        if exist('Xt','var'); Ft = Ft./100; end
    end
    vstype = 'sct';  % scatter plots
elseif strcmp(file(end-3:end),'.csv') || ...  % load data from text file
       strcmp(file(end-3:end),'.txt') || ...
       strcmp(file(end-3:end),'.dat')
    DATA   = readtable(file);
    PRJCT  = DATA.Properties.VariableNames(1    ); PRJCT = PRJCT{:};
    VNAMES = DATA.Properties.VariableNames(2:end);
    SNAMES = DATA{:    ,1};
    X      = DATA{:,2:end};
    if any(sum(X,2) > 2)  % detect that chemical data is in wt %, convert to wt fraction
        X = X./100;
        if exist('Ft','var'); Ft = Ft./100; end
        if exist('Xt','var'); Ft = Ft./100; end
    end
    vstype = 'sct';  % scatter plots
elseif strcmp(file(end-3:end),'.tif') || ...  % load image data from tiff file
       strcmp(file(end-3:end),'tiff')
    T   = Tiff(file,'r');
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
    [DGN.mx,DGN.my,n] = size(IMG);
    m   = DGN.mx*DGN.my;
    X   = reshape(IMG,m,size(IMG,3));
    X   = double(X);
    VNAMES = cell(1,n); for i=1:n; VNAMES(i) = {['BND ',int2str(i)]}; end
    SNAMES = {};
    vstype = 'img';  % image plots
else % assume data is provided as Matlab file
    load(file);
    if exist('X','var')
        if any(sum(X,2) > 2)  % detect that chemical data is in wt %, convert to wt fraction
            X = X./100;
            if exist('Ft','var'); Ft = Ft./100; end
            if exist('Xt','var'); Xt = Xt./100; end
        end
        vstype = 'sct';  % scatter plots
    elseif exist('IMG','var')
        if any(size(IMG(:)) > 1e5)  % coarsen image to avoid working with oversized data
            dft = 1;
            crs = input(['\nImage has dimensions of ',int2str(size(IMG,1)),' x ',int2str(size(IMG,2)),' x ' ,int2str(size(IMG,3)),'; apply coarsening factor (coarsen>1, dft=1):\n']);
            if isempty(crs); crs = dft; end
            IMG = IMG(1:crs:end,1:crs:end,:);
            for i = 1:length(size(IMG,3))
                IMG(:,:,i) = imnlmfilt(IMG(:,:,i));
            end
        end
        [DGN.mx,DGN.my,n] = size(IMG);
        m   = DGN.mx*DGN.my;
        X   = reshape(IMG,m,size(IMG,3));
        X   = double(X);
        VNAMES = cell(1,n); for i=1:n; VNAMES(i) = {['BND ',int2str(i)]}; end
        SNAMES = {};
        vstype = 'img';  % image plots
    end
end
if ~exist('PRJCT','var'); filesplt = split(file,'.'); PRJCT = filesplt{1}; end
X0 = X;  % store original, unprocessed data

% store number of samples and variables in diagnostics structure
[m,n]  = size(X);
DGN.m  = m; DGN.n = n; DGN.p = n-1;
DGN.Ii = (1:m).';
DGN.Ir = [];
DGN.fn = 1;

% plot unprocessed data
if exist('Ft','var')  % if true EM known
    DGN.fh(DGN.fn) = visualise(DGN.fn,{X,Xt,Ft},{'true data','noisy data','true EMs'},'Unprocessed Data',DGN,VNAMES,vstype); DGN.fn = DGN.fn+1; 
else
    DGN.fh(DGN.fn) = visualise(DGN.fn,{X},{'data'},'Unprocessed Data',DGN,VNAMES,vstype); DGN.fn = DGN.fn+1; 
    if strcmp(vstype,'img')
        DGN.fh(DGN.fn) = visualise(DGN.fn,{X},{'data'},'Unprocessed Data',DGN,VNAMES,'rgb'); DGN.fn = DGN.fn+1; 
    end
end


%% *****  REDUCE DIMENSIONALITY  ******************************************

fprintf(1,'\n*****  ANALYSE DATA \n\n')

[X,DGN]     = analyse(X,DGN);       % perform principal component analysis

fprintf(1,'\n*****  REDUCE DATA DIMENSIONS \n\n')

[Xp,Fp,Ap,DGN] = reduce(X,DGN,VNAMES); % reduce data dimensionality

PCNAMES     = cell(1,DGN.p); for i=1:DGN.p; PCNAMES(i) = {['PC ',int2str(i)]}; end % create labels for principal components

% plot reduced data
if exist('Ft','var')  % if true EM known
    DGN.fh(DGN.fn) = visualise(DGN.fn,{Xp,Xt,Ft},{'true data','reduced data','true EMs'},'Reduced Data',DGN,VNAMES,vstype); DGN.fn = DGN.fn+1; 
else
    DGN.fh(DGN.fn) = visualise(DGN.fn,{Xp},{'data'},'reduced Data',DGN,VNAMES,vstype); DGN.fn = DGN.fn+1; 
    if strcmp(vstype,'img')
        DGN.fh(DGN.fn) = visualise(DGN.fn,{Xp},{'data'},'Reduced Data',DGN,VNAMES,'rgb'); DGN.fn = DGN.fn+1; 
    end
end


%% *****  MAXIMISE INTERNAL SIMPLEX  **************************************

fprintf(1,'\n*****  FIND INTERNAL END-MEMBERS \n\n')

dft = 0;
iem = input('\n->  Do you wish to find internal end-members? yes=1, no=0 (dft) \n');
if isempty(iem); iem = dft; end

if iem
    % find internal EMs that span maximum inclusive volume in RPC space
    % based on VCA algorithm by [REF]
    [Fi,Fip,Ai,DGN] = maximise(Ap,DGN);
    
    DGN.fh(DGN.fn) = visualise(DGN.fn,{Ap(DGN.Ii,:),Fip},{'proj. data','initial EMs'},['Fitted data with ',num2str(DGN.p),' internal EMs;'],DGN,PCNAMES); DGN.fn = DGN.fn+1;
else
    Fi = ones(DGN.m,DGN.p)./DGN.p;
end


%% *****  MINIMISE EXTERNAL SIMPLEX  **************************************

fprintf(1,'\n*****  FIND EXTERNAL END-MEMBERS \n\n')

dft = 0;
eem = input('\n->  Do you wish to find external end-members? yes=1, no=0 (dft) \n');
if isempty(eem); eem = dft; end

if eem
    % find external EMs that span minimum exclusive volume in RPC space
    % based on the MINVEST algorithm by [REF]
    [Fe,Fep,Ae,DGN] = minimise(Ap,Fi,DGN);
    
    DGN.fh(DGN.fn) = visualise(DGN.fn,{Ap(DGN.Ii,:),Fep},{'proj. data','initial EMs'},['Fitted data with ',num2str(DGN.p),' external EMs;'],DGN,PCNAMES); DGN.fn = DGN.fn+1;
end


%% *****  EXAMINE RESULTS  ************************************************

fprintf(1,'\n*****  EXAMINE OUTPUT \n\n')

% report final data fit and EM compositions
if exist('Ft','var')  % if true EM known
    if iem; DGN.fh(DGN.fn) = visualise(DGN.fn,{X,Xp,Fi,Ft},{'orig. data','proj. data','internal EMs','true EMs'},['Fitted model with ',num2str(DGN.p),' internal EMs;'],DGN,VNAMES);  DGN.fn = DGN.fn+1; end
    if eem; DGN.fh(DGN.fn) = visualise(DGN.fn,{X,Xp,Fe,Ft},{'orig. data','proj. data','external EMs','true EMs'},['Fitted model with ',num2str(DGN.p),' external EMs;'],DGN,VNAMES);  DGN.fn = DGN.fn+1; end
else
    if iem; DGN.fh(DGN.fn) = visualise(DGN.fn,{X,Xp,Fi},{'orig. data','proj. data','internal EMs'},['Fitted model with ',num2str(DGN.p),' internal EMs;'],DGN,VNAMES);  DGN.fn = DGN.fn+1; end
    if eem; DGN.fh(DGN.fn) = visualise(DGN.fn,{X,Xp,Fe},{'orig. data','proj. data','external EMs'},['Fitted model with ',num2str(DGN.p),' external EMs;'],DGN,VNAMES);  DGN.fn = DGN.fn+1; end
end

if strcmp(vstype,'img')
    if iem; DGN.fh(DGN.fn) = visualise(DGN.fn,{Ai},{},['Map of ',num2str(DGN.p),' internal EM proportions;'],DGN,VNAMES,'rgb');  DGN.fn = DGN.fn+1; end
    if eem; DGN.fh(DGN.fn) = visualise(DGN.fn,{Ae},{},['Map of ',num2str(DGN.p),' external EM proportions;'],DGN,VNAMES,'rgb');  DGN.fn = DGN.fn+1; end
end


%% *****  CLUSTERING ANALYSIS  ********************************************

fprintf(1,'\n*****  ANALYSE CLUSTERING \n\n')

dft = 0;
clu = input('\n->  Do you wish to perform clustering analysis? yes=1, no=0 (dft) \n');
if isempty(clu); clu = dft; end

if clu
    % select data representation for clustering analysis
    dft = 'PC';
    cdt = input(['->  Select data representation to run clustering analysis on: \n' ...
        '    data in principal component space: PC (dft) \n' ...
        '    data in internal end-member space: iEM \n' ...
        '    data in external end-member space: eEM \n'],'s');
    if isempty(cdt); cdt = dft; end
    DGN.cluster_data = cdt;
    
    switch cdt
        case 'PC'
            Ac = Ap;
            [Ic,Fcc,DGN] = clustering(Ac,DGN);
            Fc = Fcc*DGN.PC(1:DGN.p-1,:) + DGN.meanX;
        case 'iEM'
            if iem
                Ac = Ai;
                [Ic,Fcc,DGN] = clustering(Ac,DGN);
                Fc = Fcc*Fi;
            else
                disp('no internal EM available')
            end
        case 'eEM'
            if eem
                Ac = Ae;
                [Ic,Fcc,DGN] = clustering(Ac,DGN);
                Fc = Fcc*Fe;
            else
                disp('no external EM available')
            end
    end
    
    % prepare cluster visualisation according to membership
    DATA = cell(1,DGN.c+1);
    LGND = cell(1,DGN.c+1);
    for c = 1:DGN.c
        DATA(c) = {Ac(DGN.Ii(Ic==c),:)};
        LGND(c) = {['cluster ',int2str(c)]};
    end
    DATA(end) = {Fcc};
    LGND(end) = {'cluster centroids'};
    CCNAMES   = cell(1,DGN.p); for i=1:DGN.p; CCNAMES(i) = {['CC ',int2str(i)]}; end % create labels for principal components
    
    DGN.fh(DGN.fn) = visualise(DGN.fn,DATA,LGND,['Fitted data with ',num2str(DGN.c),' cluster centroids;'],DGN,CCNAMES); DGN.fn = DGN.fn+1;
    if strcmp(vstype,'img')
        DGN.fh(DGN.fn) = visualise(DGN.fn,{Ic},{},['Map of ',num2str(DGN.c),' clusters;'],DGN,VNAMES,'img');  DGN.fn = DGN.fn+1;
    end
end


%% *****  REPORT TRAINED MODEL  *******************************************

fprintf(1,'\n*****  REPORT TRAINED MODEL \n\n')

DGN.CD = 1 - std(Xp-X).^2./std(X).^2;
disp(' ');
fprintf(1,'                         '); 
for nn=1:n
    fprintf(1,'%s   ',VNAMES{nn});
    for ll = 1:5-length(VNAMES{nn}); fprintf(1,' '); end
end; fprintf(1,'\n');
disp([    '==> DATA FIT: VAR. CD = ',num2str(DGN.CD                 ,'%1.3f   ')]);  
disp([    '              NORM CD = ',num2str(norm(DGN.CD,2)./sqrt(n),'%1.3f   ')]);
if iem
    disp(' ');
    disp(     '==> INTERNAL EM COMPOSITIONS:');
    disp(' ');
    for nn=1:DGN.n
        fprintf('%s  ',VNAMES{nn});
        for ll = 1:8-length(VNAMES{nn}); fprintf(1,' '); end
        for pp=1:DGN.p; if Fi(pp,nn)>0; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fi(pp,nn)); end
        fprintf(1,'\n');
    end
end
if eem
    disp(' ');
    disp(     '==> EXTERNAL EM COMPOSITIONS:');
    disp(' ');
    for nn=1:DGN.n
        fprintf('%s  ',VNAMES{nn});
        for ll = 1:8-length(VNAMES{nn}); fprintf(1,' '); end
        for pp=1:DGN.p; if Fe(pp,nn)>0; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fe(pp,nn)); end
        fprintf(1,'\n');
    end
end
if clu
    disp(' ');
    disp(     '==> CLUSTER CTR COMPOSITIONS:');
    disp(' ');
    for nn=1:DGN.n
        fprintf('%s  ',VNAMES{nn});
        for ll = 1:8-length(VNAMES{nn}); fprintf(1,' '); end
        for cc=1:DGN.c; if Fc(cc,nn)>0; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fc(cc,nn)); end
        fprintf(1,'\n');
    end
end
disp(' ');


%% *****  APPLY TRAINED MODEL  ********************************************

fprintf(1,'\n*****  APPLY TRAINED MODEL \n\n')

dft = 0;
app = input('\n->  Do you wish to apply trained model to other data? yes=1, no=0 (dft) \n');
if isempty(app); app = dft; end

if app;  apply;  end


%% *****  SAVE RESULTS  ***************************************************

fprintf(1,'\n*****  SAVE RESULTS \n\n')

% save output and figures
dft = 0;
sdf = input('\n->  Do you wish to save data and figures? yes=1, no=0 (dft) \n');
if isempty(sdf); sdf = dft; end

if sdf
    nowstr  = datestr(datetime('now')); nowstr(nowstr==' ') = '_';  nowstr(nowstr=='-') = []; nowstr(nowstr==':') = [];
    sdfpath = ['./out/',PRJCT,'_',nowstr,'/'];
    if ~exist(sdfpath,'dir'); mkdir(sdfpath); end
    save([sdfpath,'out']);
    for fn = 1:length(DGN.fh)
        if ishandle(DGN.fh(fn))
            set(DGN.fh(fn),'PaperSize',[24 18]);
            print(DGN.fh(fn),[sdfpath,'fig_',int2str(fn)],'-dpdf','-r300','-loose');
            savefig(DGN.fh(fn),[sdfpath,'fig_',int2str(fn)]);
        end
    end
    if exist('DGNa','var')
        for fn = 1:length(DGNa.fh)
            if ishandle(DGNa.fh(fn))
                set(DGNa.fh(fn),'PaperSize',[24 18]);
                print(DGNa.fh(fn),[sdfpath,'fig_',int2str(fn)],'-dpdf','-r300','-loose');
                savefig(DGNa.fh(fn),[sdfpath,'fig_',int2str(fn)]);
            end
        end
    end
end
