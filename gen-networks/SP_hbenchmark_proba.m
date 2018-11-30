%
% function [p,S] = SP_hbenchmark_proba(COM, rho, kbar)
%
% Computes the probabilities of link existence necessary to 
% create Sales Pardo graphs
%
% Inputs : 
% - COM, the multiscale structure of size (N)x(number of
% structural levels). For a structural level i, COM(j,i) is
% the number of the partition that contains node j  
% - rho and kbar, two parameters of Sales Pardo graphs.
%
% Outputs :
% - p, the probabilities of link existence 
% corresponding to the multiscale community structure coded 
% in COM and parameters rho and kbar
% - S corresponds to the number of nodes in each level of scale.
%
% This function follows the definition of a multiscale benchmark graph one
% can find in the supplementary material of "Extracting the hierarchical 
% organization of complex systems", Sales-Pardo and Guimera and Moreira and 
% Amaral, PNAS (2007). Nevertheless, these authors are not 
% responsible for any error in this implementation. Errors should 
% be addressed to nicolas DOT tremblay AT ens-lyon DOT fr

% This file is part of the MSCD_Wav toolbox (MultiScale 
% Community Detection using Wavelets toolbox)
% Copyright (C) 2014, Nicolas Tremblay. 
%
% The MSCD_Wav toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The MSCD_Wav toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the MSCD_Wav toolbox.  If not, see <http://www.gnu.org/licenses/>.

function [p,S] = SP_hbenchmark_proba(COM, rho, kbar)

% small verif:
for l=1:size(COM,2)
    [x,y]=hist(COM(:,l),max(COM(:,l)));
    if std(x)~=0
        error('SP_hbenchmark.m: problem in COM: it cannot code for a SP graph')
    else
        S(l)=x(1);
    end
end

S=[S(1)*(max(COM(:,1))-1), S];
for l=2:size(COM,2)
    S(l)=S(l)-S(l+1);
end
S(end)=S(end)-1;

l=1;
p(1)=((rho^size(COM,2))/((1+rho)^size(COM,2)))*(kbar/S(1));

for l=2:size(COM,2)+1
    p(l)=((rho^(size(COM,2)-(l-1)))/((1+rho)^(size(COM,2)-(l-1)+1)))*(kbar/(S(l)));
end
