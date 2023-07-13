function [w, T] = Sparse_Diff(G)
%UNTITLED2 Summary of this function goes here
%  input:
%  G--the Adjacent matrix of a graph

n=size(G,1);
w=(rand(n^2,1)+0.7)/2;
B=zeros(n,n);

k=1;
for i=1:n
    for j=i+1:n
        if G(i,j)>0
            B(i,j)=w(k);
            k=k+1;
        end
    end
end

B=B+B';

%create the weight degree matrix
D=diag(sum(B,2));


T=pinv(sqrt(D))*B*pinv(sqrt(D));
w=w(1:k-1);
end





