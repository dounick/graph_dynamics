function [y] = Fun4(a,I,S,t,T)
% construct target minimization function 
%a: weight variables, w: original weight, I: dirac sources locations
%J:sampling locations, t: sampling intervals T: diffusion matrix



%  Diffusion matrix

W=Diff_Matrix(a,T);


% initialization of target function
y=0;

% sampling trajectories 

for i=1:length(t)
  l=t(i);
  B=T^l-W^l;
  y=y+norm(B(S,I));
  
  %+(2*rand(n,m)-1)*10^(-10));
end


end