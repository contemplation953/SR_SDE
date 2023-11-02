clc;clear;


data_path=sprintf('%s/5-dof/data',pwd);
ans_path=sprintf('%s/Theoretical_Value/data',pwd);

data_index=1;
data=load(sprintf('%s/H_%d.mat',data_path,data_index));
pst=0.01:0.01:25;
pdf_bound=1e-6;
monment_count=data.monment_count;
idnex=find(monment_count/sum(monment_count)>=pdf_bound);
begin_idnex=idnex(1);
end_index=idnex(end);
pst=pst(begin_idnex:end_index);
moment1=data.moment1(begin_idnex:end_index);
moment2=data.moment2(begin_idnex:end_index);

pst1=0.001:0.001:30;
data1=load(sprintf('%s/analyse_%d.mat',ans_path,data_index));
M_H=data1.M_H;
Sigma_H_square=data1.Sigma_H_square;


data2=load(sprintf('%s/global_params.mat',data_path));
c=data2.x;

D_list=0.005:0.0005:0.035;
x=D_list(data_index);


m_x=@(y) c(1)+c(2).*x+c(3).*y+c(4).*x.*y+c(5).*y.^2;
s_x=@(y) c(6)+c(7).*x+c(8).*y+c(9).*x.^2+c(10).*x.*y+c(11).*y.^2+c(12).*x.^2.*y+c(13).*x.*y.^2+c(14).*y.^3;
pst2=0.5:0.5:20;
sim_m=m_x(pst2);
sim_s=(s_x(pst2)).^2;


figure(1);
plot(pst,moment1,'kO',pst1,M_H,'r',pst2,sim_m,'b+');
legend('1d-moment','analyse','sim_res');
xlabel('H');
ylabel('MSE');
grid on;
title('m(H)');

figure(2);
plot(pst,moment2,'kO',pst1,Sigma_H_square,'r',pst2,sim_s,'b+');
legend('1d-moment','analyse','sim_res');
grid on;
xlabel('H');
ylabel('MSE');
title('\sigma^2(H)');

