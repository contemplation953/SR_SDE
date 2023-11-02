clc;clear;

data_path=sprintf('%s/N-DOF/10_dof_data',pwd);
data_index=33;
D_list=0.005:0.0005:0.035;
data=load(sprintf('%s/H_%d.mat',data_path,data_index));
pst=0.01:0.01:30;
pdf_bound=1e-2;
p_data=data.p_data;
idnex=find(p_data>pdf_bound);
begin_idnex=idnex(1);
end_index=idnex(end);
pst=pst(begin_idnex:end_index);
moment1=data.moment1(begin_idnex:end_index);
moment2=data.moment2(begin_idnex:end_index);


data2=load(sprintf('%s/global_params.mat',data_path));
c=data2.x;
pst2=0.5:0.5:20;

d=D_list(data_index);
sim_m_fun=@(h) c(1)+c(2)*d+c(3)*h+c(4)*d.^2+c(5)*d.*h+c(6)*h.^2;
sim_s_fun=@(h) (c(7)+c(8)*d+c(9)*h+c(10)*d.^2+c(11)*d.*h+c(12)*h.^2+c(13)*d.^2.*h+c(14)*d.*h.^2+c(15)*h.^3).^2;

sim_m=sim_m_fun(pst2);
sim_s=sim_s_fun(pst2);

figure(1);
plot(pst,moment1,'kO',pst2,sim_m,'b+');
legend('1d-moment','sim_res');
title('M(H)');

figure(2);
plot(pst,moment2,'kO',pst2,sim_s,'b+');
legend('2d-moment','sim_res');
title('\sigma^2(H)');

