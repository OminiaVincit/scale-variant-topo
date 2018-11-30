%
% function A = SP_hbenchmark(COM, rho, kbar)
%
% Creates a realisaton of a Sales-Pardo graph.
%
% Inputs : 
% - COM, the multiscale structure of size (N)x(number of
% structural levels). For a structural level i, COM(j,i) is
% the number of the partition that contains node j  
% - rho and kbar, two parameters of Sales Pardo graphs.
%
% Output :
% - the adjacency matrix of the graph
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

function A = SP_hbenchmark(COM, rho, kbar)

N=size(COM,1);
[p,S] = SP_hbenchmark_proba(COM, rho, kbar);

P=zeros(N,N);

for l=1:size(COM,2)
    S=repmat(COM(:,l),1,max(COM(:,l)))==repmat([1:max(COM(:,l))],N,1);
    for k=1:max(COM(:,l))
        P(S(:,k),S(:,k))=p(l+1);
    end
end

P(P==0)=p(1);
rng('shuffle');
A=ceil(P-rand(N));
A=triu(A,1);A=A+A';
