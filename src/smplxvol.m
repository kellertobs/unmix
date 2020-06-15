% smplxvol:  unmix subroutine to calculate proxy simplex volume 
%
%   [f] = smplxvol(F)
%
%   The routine calculates a proxy to track the volume of a (p-1)-dimensional 
%   simplex used in the objective function of the optimisation for
%   identifying external mixing end-members. The simplex volume proxy is
%   calculated as
%
%   V = sqrt(abs(det(F*F.')))
%
%   for the non-square end-member matrix F [p x p-1] in RCP space.
%
%   F      : input end-member compositions in RCP space
%
%   V      : output volume proxy for (p-1)-dimensional simplex
%
% created  : 2020-05-05  Tobias Keller, University of Glasgow
% license  : GNU General Public License v3.0


function V = smplxvol(F)

F = [F,ones(size(F,1),1)];
V = sqrt(abs(det(F.'*F)));

end  % end function
