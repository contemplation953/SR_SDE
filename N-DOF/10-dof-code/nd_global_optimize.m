clc;clear

data_path=sprintf('%s/N-DOF/10_dof_data',pwd);

data=load(sprintf('%s/10dofglobal_init_params.mat',data_path));
p_m=data.params_m;
p_s=data.params_s;
H_list=0.01:0.01:30;
D_list=0.005:0.0005:0.035;
pdf_matrix=zeros(61,length(H_list));
moment1_matrix=zeros(61,length(H_list));
moment2_matrix=zeros(61,length(H_list));
coeff=[p_m,p_s];
    
tic;
for i=1:61
    data=load(sprintf('%s/H_%d.mat',data_path,i));
    pdf_matrix(i,:)=data.p_data(1:length(H_list));
    moment1_matrix(i,:)=data.moment1(1:length(H_list));
    moment2_matrix(i,:)=data.moment2(1:length(H_list));
end

%options = optimset('PlotFcns',@optimplotfval,'MaxFunEvals',1e9,'MaxIter',1e2,'TolFun ',1e-4);
options = optimset('MaxFunEvals',1e9,'MaxIter',2e3,'TolFun ',1e-4);
fun=@(coeff) nd_global_cost(coeff,pdf_matrix,H_list,D_list,moment1_matrix,moment2_matrix);
%一次大概6s
[x,val]=fminsearch(fun,coeff,options);


toc;   
save(sprintf('%s/global_params.mat',data_path),'x');

