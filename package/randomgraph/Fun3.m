function [y] = Fun3(a,w,I,S,t,L)
% construct target minimization function 
%a: weight variables, w: original weight, I: dirac sources locations
%J:sampling locations, t: sampling intervals



%  Diffusion matrix
W=Reg_Diff_Matrix(w,L);
A=Reg_Diff_Matrix(a,L);



% initialization of target function
y=0;

% sampling trajectories 

for i=1:length(t)
  l=t(i);
  B=A^l-W^l;
  y=y+norm(B(S,I));
  
  %+(2*rand(n,m)-1)*10^(-10));
end


end