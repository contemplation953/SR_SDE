clc;clear;

index=1;
data_path=sprintf('%s/5-dof/data',pwd);
ans_path=sprintf('%s/Theoretical_Value/data',pwd);

h=0.01;
H_list=h:h:25-h;

data=load(sprintf('%s/H_%d.mat',data_path,index));
moment1=data.moment1(1:length(H_list));
moment2=sqrt(data.moment2(1:length(H_list))); 
monment_count=data.monment_count;
w=monment_count/sum(monment_count);

data1=load(sprintf('%s/analyse_%d.mat',ans_path,index));
ans_h_list=0.001:0.001:30;
M_H=data1.M_H;
Sigma_H_square=sqrt(data1.Sigma_H_square);

p1=polyfit(H_list,moment1(1:length(H_list)),3);
p2=polyfit(H_list,moment2(1:length(H_list)),5);

pst=0.1:0.1:20;
y1=polyval(p1,pst);
y2=polyval(p2,pst);


figure(1);
plot(ans_h_list,M_H,'b-',H_list,moment1(1:length(H_list)),'k.',pst,y1,'r-');

figure(2);
plot(ans_h_list,Sigma_H_square,'b-',H_list,moment2(1:length(H_list)),'k.',pst,y2,'r-');


