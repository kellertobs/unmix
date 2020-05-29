% minimise:  unmix subroutine to minimise volume of data-external simplex
%
%   [F,Fp,A] = minimise(X0,A0,F0,DGN)
%
%   The routine calculates external mixing end-members and corresponding 
%   mixing abundances that best fit the data as projected into RPC space by
%   minimising the volume of the data-external simplex spanned by the end-
%   members while maintaining non-negative abundances and end-member 
%   compositions to a set tolerance. The routine is based on the MINVEST 
%   algorithm developed by [REF] and returns end-member coordinates and 
%   abundances in RPC and MFC space. 
%
%   A      : input mixing abundances to fit data in RCP space
%   F      : input starting guess for end-member compositions in FCM space
%   DGN    : input structure containing various diagnostics
%
%   F      : output calculated end-member compositions in FCM space
%   Fp     : output calculated end-member compositions in RCP space
%   A      : output corresponding mixing abundances in FCM space
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function [F,Fp,A,DGN] = minimise(A,F,DGN)

p = DGN.p;  % #end-members

% set tolerances for EM selections
dft = [0.01,0.01];
tol = input(['->  Adjust endmember fitting tolerances as list [Atol,Ftol] \n' ...
             '    Atol: tolerance for negative abundances    (dft = ',num2str(dft(1)),') \n' ...
             '    Ftol: tolerance for negative EM components (dft = ',num2str(dft(2)),') \n']);
if isempty(tol); tol = dft; end
Atol   = tol(1);
Ftol   = tol(2);

Fp     = (F - DGN.meanX)/DGN.PC(1:p-1,:);

opts   = optimoptions('fmincon','Display','iter');
[Fp,~] = fmincon(@(Fp)log(smplxvol(Fp)),Fp,[],[],[],[],[],[],@(Fp)nonneg(Fp,A(DGN.Ii,:),Atol,Ftol,DGN),opts);

F  = Fp*DGN.PC(1:p-1,:) + DGN.meanX;  % get min-vol. external EMs in FCM space

A  = [A,ones(size(A,1),1)]/[Fp,ones(size(Fp,1),1)];  % get mixing abundances in FCM space

DGN.minV = smplxvol(Fp);  % get minimised data-exclusive simplex volume

end  % end function

