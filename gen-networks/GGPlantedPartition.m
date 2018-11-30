function [A,V0]=GGPlantedPartition(NN,pi,pe,Diag)
% function [A,V0] = GGPlantedPartition(NN,pi,pe,Diag)
% Generation of a planted partition graph
%
% Creates a planted partition graph, returns its adjacency matrix.
% This is the classical planted partition graph, with unequal-sized 
% partitions. For detals see A.Condon and R.M. Karp. "Algorithms 
% for graph partitioning on the planted partition model.",  
% Randomization, Approximation, and Combinatorial Optimization. 
% Algorithms and Techniques. Springer Berlin Heidelberg, 1999. 
%
% INPUT:
% NN:    vector of  community boundaries
% pi:    internal edge probability
% pe:    external edge probability
% Diag:  if Diag=1, use self-loops; if Diag=0, don't use self-loops
%
% OUTPUT:
% A  adajcency matrix (N-by-N)
% V0 classification vector (N-by-1)
% W  permutation matrix (N-by-N) to put the nodes in order
%
%EXAMPLE
% [A,V0]=GGPlantedPartition([0 10 20 30 40],0.9,0.1,0);
%
K=length(NN)-1;
N=NN(K+1);

A0=eye(N);
for k=1:K
    N1=NN(k)+1;
    N2=NN(k+1);
    A0(N1:N2,N1:N2)=1;
	N0=N2-N1;
	MiMax(k)=N0*(N0-1)/2;
end

A=eye(N);
for n1=1:N
    for n2=n1+1:N
        if A0(n1,n2)==1 & rand(1)<pi; A(n1,n2)=1; A(n2,n1)=1; end
        if A0(n1,n2)==0 & rand(1)<pe; A(n1,n2)=1; A(n2,n1)=1; end
    end
end
%A=A-eye(N);

i=1;
for k=1:K
	V0(NN(k)+1:NN(k+1),1)=k;
	MI(k,1)=sum(sum(A(NN(k)+1:NN(k+1),NN(k)+1:NN(k+1))))/2;
	MT(k,1)=sum(sum(A(NN(k)+1:NN(k+1),:)))/2;
	ME(k,1)=MT(k,1)-MI(k,1);
	aa(k,1)=ME(k)/MI(k);
	bb(k,1)=MI(k)/MiMax(k);
end

if Diag~=0
	for n=1:N; A(n,n)=1; end
else	
	for n=1:N; A(n,n)=0; end
end
