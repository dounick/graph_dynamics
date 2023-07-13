function [T] = Cir_Diff_Matrix(w)
% construct diffusion matrix for circular graph
% If G is a circular graph, we associate each edge a weight and w
% is the weight vector. 

% number of edges
n=length(w);

%initialization of diffusion matrix
T=zeros(n,n);


T(1,2)=w(1);
T(1,n)=w(n);
for i=2:n-1
    T(i,i+1)=w(i);
end

% Weight matrix 
T=T+T';
% Weight degree matrix
D=diag(sum(T,2));

% Diffusion matrix
T=pinv(sqrt(D))*T*pinv(sqrt(D));
end

