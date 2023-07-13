A=[];
b=[];
Aeq=[];
beq=[];

nv=64; Kreg=2; p=0;
[G]=erdosRenyi(nv,p,Kreg);
[w, T] = erdosRenyi_diffusion(G);
% creat a weighted graph using Erdos-Renyi graph model
%nv=64;p=0.05;
%rand('seed',100);
%G = rand(nv,nv) < p;
%G = triu(G,1);
%G = G + G';
%[w,T]=Sparse_Diff(G);

%Time=8:16;
index1=randperm(nv);
I=index1(floor(0.2*nv));
Exp_times=1;
err=zeros(Exp_times);
  
  
   for k=1:Exp_times
        % choose source locations
        J=1:2:nv;
        fun=@(x)Fun4(x,J,I,64,T);
        options = optimset('MaxIter',3000000,'MaxFunEvals',3000000,'Display','iter');
        x=fmincon(fun,rand(length(w),1),A,b,Aeq,beq,zeros(length(w),1)*0.000001,ones(length(w),1)*0.999999,[],options);
        err(k)=norm(T-Diff_Matrix(x,T));
   end

   
 
   

 

 
   


 