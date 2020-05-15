% nonneg:  unmix subroutine to impose non-negativity constraints on optimisation
%
%   [c,ceq] = nonneg(F,A,Atol,Ftol,minA,maxA,PC,mX0)
%
%   The routine calculates end-member compositions and abundances in
%   appropriate data space and constrains them to non-negative values or 
%   up to the set tolerances -Ftol, and -Atol, respectively.
%
%   F      : input current guess of external end-member compositions in RPC
%   A      : input 
%   Atol   : input tolerance of non-negative end-member abundances set by user
%   Ftol   : input tolerance of non-negative end-member composition set by user
%   PC     : input principal component compositions in FMC space
%   mX0    : input variable-wise mean of sample compositions
%
%   c      : output inequality used to constrain optimisation (c <= 0)
%   ceq    : output equality used to constrain optimisation (ceq <= 0)
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function [c,ceq] = nonneg(F,A,Atol,Ftol,DGN)

Fp  = [F,ones(size(F,1),1)];
Ap  = [A,ones(size(A,1),1)];
A   = Ap/Fp;

F   = F*DGN.PC(1:DGN.p-1,:) + DGN.meanX;

c   = [-(A(:) + Atol); -(F(:) + Ftol)];
ceq = [];

end  % end function
