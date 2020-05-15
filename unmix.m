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
clear variables;       % clear workspace
warning('off','all');  % turn off warning
addpath ./src ./data   % add paths to source code and data directories


%% *****  LOAD DATA  ******************************************************

% set data file to read
dft  = 'DATA';
file = input(['->  What dataset would you like to process? \n', ...
              '    Give full name of datafile (dft = DATA): \n'],'s');
if isempty(file); file = dft; end

% read file to load dataset
if strcmp(file(end-3:end),'.csv') || ...  % load data from text file
   strcmp(file(end-3:end),'.txt') || ...
   strcmp(file(end-3:end),'.dat')
    DATA   = readtable(file);
    PRJCT  = DATA.Properties.VariableNames(1    ); PRJCT = PRJCT{:};
    VNAMES = DATA.Properties.VariableNames(2:end);
    SNAMES = DATA{:    ,1};
    X      = DATA{:,2:end};
else % assume data is provided as Matlab file
    load(file);
end

% store number of samples and variables in diagnostics structure
[m,n]  = size(X);
DGN.m  = m; DGN.n = n; DGN.p = n-1;
DGN.Ii = (1:m).';
DGN.Ir = [];

% plot unprocessed data
if exist('Ft','var')  % if true EM known
    f1 = visualise(1,{X,Xt,Ft},{'true data','noisy data','true EMs'},'Unprocessed Data',DGN,VNAMES);
else
    f1 = visualise(1,{X},{'data'},'Unprocessed Data',DGN,VNAMES);
end
pause


%% *****  REDUCE DIMENSIONALITY  ******************************************

[X,DGN]     = analyse(X,DGN);       % perform principal component analysis

[Xp,Ap,DGN] = reduce(X,DGN,VNAMES); % reduce data dimensionality


%% *****  CLUSTERING ANALYSIS  ********************************************

[Ic,Fc,Fcp,DGN] = clustering(Ap,DGN);

% prepare cluster visualisation according to membership
DATA = cell(1,DGN.p+1);
DATX = cell(1,DGN.p+1);
LGND = cell(1,DGN.p+1);
for c = 1:DGN.p
    DATA(c) = {Ap(DGN.Ii(Ic==c),:)};
    DATX(c) = {Xp(DGN.Ii(Ic==c),:)};
    LGND(c) = {['cluster ',int2str(c)]};
end
DATA(end) = {Fcp};
DATX(end) = {Fc};
LGND(end) = {'cluster centroids'};
PCNAMES = cell(1,DGN.p); for i=1:DGN.p; PCNAMES(i) = {['PC',int2str(i)]}; end

f2 = visualise(2,DATA,LGND,['Fitted data with ',num2str(DGN.p),' cluster centroids;'],DGN,PCNAMES);
pause


%% *****  MAXIMISE INTERNAL SIMPLEX  **************************************

% find internal EMs that span maximum inclusive volume in RPC space
% based on VCA algorithm by [REF]
[Fi,Fip,Ai,DGN] = maximise(Ap,Fcp,DGN);

f3 = visualise(3,{Ap(DGN.Ii,:),Fip},{'proj. data','initial EMs'},['Fitted data with ',num2str(DGN.p),' internal EMs;'],DGN,PCNAMES);
pause


%% *****  MINIMISE EXTERNAL SIMPLEX  **************************************

% find external EMs that span minimum exclusive volume in RPC space
% based on the MINVEST algorithm by [REF]
[Fe,Fep,Ae,DGN] = minimise(Ap,Fi,DGN);

f4 = visualise(4,{Ap(DGN.Ii,:),Fep},{'proj. data','initial EMs'},['Fitted data with ',num2str(DGN.p),' external EMs;'],DGN,PCNAMES);
pause


%% *****  REPORT RESULTS  *************************************************

% report final data fit and EM compositions
if exist('Ft','var')  % if true EM known
    f5 = visualise(5,DATX,LGND,['Fitted model with ',num2str(DGN.p),' clusters;'],DGN,VNAMES);
    f6 = visualise(6,{X,Xp,Fi,Ft},{'orig. data','proj. data','internal EMs','true EMs'},['Fitted model with ',num2str(DGN.p),' internal EMs;'],DGN,VNAMES);
    f7 = visualise(7,{X,Xp,Fe,Ft},{'orig. data','proj. data','external EMs','true EMs'},['Fitted model with ',num2str(DGN.p),' external EMs;'],DGN,VNAMES);
else
    f5 = visualise(5,DATX,LGND,['Fitted model with ',num2str(DGN.p),' clusters;'],DGN,VNAMES);
    f6 = visualise(6,{X,Xp,Fi},{'orig. data','proj. data','internal EMs'},['Fitted model with ',num2str(DGN.p),' internal EMs;'],DGN,VNAMES);
    f7 = visualise(7,{X,Xp,Fe},{'orig. data','proj. data','external EMs'},['Fitted model with ',num2str(DGN.p),' external EMs;'],DGN,VNAMES);
end

DGN.CD = 1 - std(Xp-X).^2./std(X).^2;
[~,i]  = sort(Fe(:,1));
disp(' ');
disp('==============================================================');
disp(' ');
fprintf(1,'                         '); 
for nn=1:n
    fprintf(1,'%s   ',VNAMES{nn});
    for ll = 1:5-length(VNAMES{nn}); fprintf(1,' '); end
end; fprintf(1,'\n');
disp([    '==> DATA FIT: VAR. CD = ',num2str(DGN.CD                 ,'%1.3f   ')]);  
disp([    '              NORM CD = ',num2str(norm(DGN.CD,2)./sqrt(n),'%1.3f   ')]);    fprintf(1,'\n');
disp(     '==> INTERNAL EM COMPOSITIONS:');
disp(' ');
for nn=1:DGN.n
    fprintf('%s  ',VNAMES{nn}); 
    for ll = 1:8-length(VNAMES{nn}); fprintf(1,' '); end
    for pp=1:DGN.p; if Fi(i(pp),nn)<10; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fi(i(pp),nn)); end
    fprintf(1,'\n'); 
end
disp(' ');
disp(     '==> EXTERNAL EM COMPOSITIONS:');
disp(' ');
for nn=1:DGN.n
    fprintf('%s  ',VNAMES{nn}); 
    for ll = 1:8-length(VNAMES{nn}); fprintf(1,' '); end
    for pp=1:DGN.p; if Fe(i(pp),nn)<10; fprintf(1,' '); end; fprintf(1,'%2.3f  ',Fe(i(pp),nn)); end
    fprintf(1,'\n'); 
end
disp(' ');