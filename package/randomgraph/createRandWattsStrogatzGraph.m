function [T,w,L] = createRandWattsStrogatzGraph(G)
%This function generate diffusion matrix of ErdosRenyi graph
% Inputs:
% G: structure matrix for graph 
% 
%  Outputs:
%  w: weights for edges 
%  T: diffusion matrix 


% give random weights for edges
w1=(rand(size(G.Edges,1),1)+0.7)/2;
L=adjacency(G);
%create symmetric adjacent matrix 
A=triu(full(L));

B=A;
k=1;
for i=1:size(A,1)
    for j=i+1:size(A,2)
        if B(i,j)>0
            B(i,j)=w1(k);
            k=k+1;
        end
    end
end

B=B+B';
w=w1(1:k-1);
%create the weight degree matrix
D=diag(sum(B,2));


T=pinv(sqrt(D))*B*pinv(sqrt(D));

end