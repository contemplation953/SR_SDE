clc;clear;

mh_p=3;
sh_p=3;
root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);
data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));

i=1;
ans_path=sprintf('%s/Theoretical_Value/data',pwd);
data1=load(sprintf('%s/analyse_%d.mat',ans_path,i));
M_H=data1.M_H;
Sigma_H_square=data1.Sigma_H_square;

H_list=0.01:0.01:25;



data2=load(sprintf('%s/H_%d.mat',data_path,i));
moment1=data2.moment1;
moment2=data2.moment2;

pst=0.001:0.001:30;
plot(pst,Sigma_H_square,'k-',H_list(1:2499),moment2(1:2499),'b+');
plot(pst,M_H,'k-',H_list(1:2499),moment1(1:2499),'b+');