% maximise:  unmix subroutine to maximise volume of data-internal simplex
%
%   [F,Fp,A,DGN] = maximise(X,A,F,DGN)
%
%   The routine calculates internal mixing end-members and corresponding 
%   mixing abundances where end-members coincide with mutually extreme data 
%   points in RPC space. This is achieved by identifying the set of outer
%   data points that maximise the volume of the data-inclusive simplex they
%   span. The routine is based on the VCA algorithm developed by [REF] and 
%   returns end-member coordinates and abundances in RPC and MFC space. 
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


function [F,Fp,A,DGN] = maximise(A,F,DGN)

p = DGN.p;

% create testing array for extracting extreme samples
At = [A(DGN.Ii,:),ones(length(DGN.Ii),1)].';

maxV = smplxvol(F);  % initialise volume for maximisation
for k=1:3*p  % perform VCA 3 x p times to ensure max-vol solution is found
    
    ind    = zeros(1,p);
    B      = zeros(p,p);
    B(p,1) = 1;
    
    % test for extreme samples starting from random number
    for i = 1:p
        w = rand(p,1);
        f = w - B*pinv(B)*w;
        f = f / sqrt(sum(f.^2));
        
        [~,ind(i)] = max(abs(f'*At));
        B(:,i)     = At(:,ind(i));
    end
    
    % check if new simplex spans maximum volume
    V = smplxvol(A(DGN.Ii(ind),:));
    if V >= maxV
        maxV = V;
        maxi = DGN.Ii(ind);
    end
end

Fp = A(maxi,:);  % extract max-vol. internal EMs in RPC space
F  = Fp*DGN.PC(1:p-1,:) + DGN.meanX;  % get min-vol. external EMs in FCM space

A  = [A,ones(size(A,1),1)]/[Fp,ones(size(Fp,1),1)];  % get mixing abundances in FCM space

DGN.maxV = maxV; % get maximised data-inclusive simplex volume
DGN.maxi = maxi; % get indices for mutually extreme samples

end  % end function
