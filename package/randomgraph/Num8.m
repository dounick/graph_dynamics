% choose only one delta_i randomly
% Choose J=V, no spatial downsampling 
A=[];
b=[];
Aeq=[];
beq=[];
m=1;
%number of nodes 
n=20;
% success ratio
ratio=zeros(1,1);
times=5;
deg=3;
L=createRandRegGraph(n,deg);

flag=0;
err=zeros(times,1);
for j=1
r=(rand(n(j)*deg/2,1)+0.7)/2;
for i=1:times
index1=randperm(n(j));
I=index1(1:floor(0.2*n(j)));  
J=1:m:n;
fun=@(x)Fun3(x,r,I,J,n(j),L);
options = optimset('MaxIter',3000000,'MaxFunEvals',3000000,'Display','iter');
x=fmincon(fun,rand(length(r),1),A,b,Aeq,beq,zeros(length(r),1)*0.0001,ones(length(r),1)*0.999999,[],options);
err(i,j)=norm(Reg_Diff_Matrix(x,L)-Reg_Diff_Matrix(r,L));

if err(i,j)<10^(-4)
    flag=flag+1;
end
end

ratio(j)=flag./times;
flag=0;



end
