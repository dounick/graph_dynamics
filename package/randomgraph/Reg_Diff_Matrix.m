function [T] = Reg_Diff_Matrix(w,A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

B=triu(full(A));
k=1;

for i=1:size(A,1)
    for j=i+1:size(A,2)
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

