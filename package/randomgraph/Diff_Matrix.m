function [T] = Diff_Matrix(w,T)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

B=triu(T);
k=1;

for i=1:size(T,1)
    for j=i+1:size(T,2)
        if B(i,j)>0
            B(i,j)=w(k);
            k=k+1;
        end
    end
end




B=B+B';
% Weight degree matrix
D=diag(sum(B,2));

% Diffusion matrix
T=pinv(sqrt(D))*B*pinv(sqrt(D));


end