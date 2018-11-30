function [A,V0]=GGReadEdgeList(EdgeFile,PartitionFile,Diag)
% function [A,V0] = GGReadEdgeList(EdgeFile,PartitionFile,Diag)
% Graph Generation from an edge list
%
% Reads a graph from ede list file and a partition from partition 
% membership files 
%
% INPUT:
% EdgeFile         filename (string) contains an M-by-2 list  
%                  of edges (as node pairs)
% PartitionFile    filename (string) contains an N-by-1 list  
%                  of cluster number to which each node belongs
% Diag:            if Diag=1/0, use/don't use self-loops
% 
% OUTPUT:
% A                adajcency matrix (N-by-N)
% V0               classification vector (N-by-1)
%
%EXAMPLE
% [A,V0]=GGReadEdgeList('e01.txt','v01.txt',0);
%
E=load(EdgeFile);
V0=load(PartitionFile);
M=size(E,1);
N=max(max(E));
if Diag==0
	A=zeros(N,N);
else
	A=eye(N);
end

for m=1:M
	u=E(m,1);
	v=E(m,2);
	A(u,v)=1;
	A(v,u)=1;
end
