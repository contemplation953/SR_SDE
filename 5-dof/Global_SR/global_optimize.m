clc;clear;close all;

data_path=sprintf('%s/5-dof/data',pwd);

h=0.01;
H_list=0.01:h:20-h;
D_list=0.005:0.0005:0.035;
pdf_matrix=zeros(61,length(H_list));

data=load(sprintf('%s/global_init_params.mat',data_path));
coeff=data.x;
mhv1=data.mhv1;
shv1=data.shv1;

tic;
for i=1:61
    data=load(sprintf('%s/H_%d.mat',data_path,i));
    pdf_matrix(i,:)=data.p_data(1:length(H_list));
end

options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e9,'MaxIter',1e4,'TolFun ',1e-5);
%options = optimset('MaxFunEvals',1e9,'MaxIter',1e5,'TolFun ',1e-5);
fun=@(coeff) global_cost(coeff,pdf_matrix,H_list,D_list,mhv1,shv1);
%一次大概6s
[x,val]=fminsearch(fun,coeff,options);
val

toc;   
save(sprintf('%s/global_params.mat',data_path),'x');

