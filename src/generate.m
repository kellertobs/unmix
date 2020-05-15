% generate:  unmix subroutine to generate artificial data for testing.
% 
%   The script randomly generates datasets for specified number of endmembers,
%   compositional variables, and samples. Each sample is found as a mixture 
%   of endmembers with a small, random 'measurment error' added. A small
%   number of samples is additionally shifted to mimick outlier data. The
%   generated data is saved to in Matlab file format suitable for loading
%   into the main rs-pva routine.
%
% created :    2020-03-21  Tobias Keller, University of Glasgow
% license :    GNU General Public License v3.0


% set parameters for generating data
m = 80;   % number of samples
n = 10;   % number of variables
p = 3;    % number of endmembers
a = 1;    % amplitude of added noise [wt %]

% generate true endmember compositions
Ft = rand(p,n); 
Ft = Ft./sum(Ft,2).*100;  % in wt %

% generate true mixing proportions
At = rand(m,p);
At = At./sum(At,2);  % in wt fract.

% generate true sample compositions
Xt = At*Ft;

% perturb sample compositions with normally distributed noise
X = Xt + a.*(rand(size(Xt))-0.5);

% add 5% outlier points
o       = floor(0.05*m);
io      = randi(m,o,1);
X(io,:) = X(io,:) .* (1-(randi(2,o,n)-1).*0.25);

% set sample, variable, and project names
SNAMES = cell(1,m); for i=1:m; SNAMES(i) = {['PC',int2str(i)]}; end
VNAMES = cell(1,m); for i=1:n; VNAMES(i) = {['PC',int2str(i)]}; end
PRJCT = ['Data generated ',datestr(now)]; 
PRJCT(PRJCT==' ' | PRJCT=='-') = '_';

DGN.m = m; DGN.n = n; DGN.p = p;
visualise({Xt,X,Ft},{'true data','noisy data','true EMs'},['Synthetic data for ',num2str(p),' EMs'],DGN,VNAMES)
        
save('data/DATA.mat','X','Xt','At','Ft','SNAMES','VNAMES','PRJCT');

% end of script
