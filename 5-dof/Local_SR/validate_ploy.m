clc;clear;

root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);
ans_path=sprintf('%s/Theoretical_Value/data',root_path);

data_index=1;
data=load(sprintf('%s/H_%d.mat',data_path,data_index));
h=0.01;
H_list=h:h:25-h;

moment1=data.moment1(1:length(H_list));
moment2=data.moment2(1:length(H_list));

pst1=0.001:0.001:30;
data1=load(sprintf('%s/analyse_%d.mat',ans_path,data_index));
M_H=data1.M_H;
Sigma_H_square=data1.Sigma_H_square;


mh_p=4;
sh_p=5;
root_path=pwd();
data2=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));
coeff=data2.x_list(data_index,:);

coeff1=coeff(1:mh_p+1);
coeff2=coeff(mh_p+2:end);
[drift,diffusion]=FuncFactory(coeff1,coeff2);

pst2=0.5:0.5:20;
sim_m=drift(pst2);
sim_s=(diffusion(pst2)).^2;


figure(1);
plot(H_list(1:500),moment1(1:500),'kO',pst1,M_H,'r',pst2,sim_m,'b+');
lgd=legend('simulation data','theoretical results','predictive results');
lgd.FontSize = 20;
xlabel('H','FontSize',20);
ylabel('MSE','FontSize',20);
grid on;
title('m(H)','FontSize',20);
saveas(gcf,sprintf('%s/res/pol_m_%d.png',root_path,mh_p));


figure(2);
plot(H_list(1:1000),moment2(1:1000),'kO',pst1,Sigma_H_square,'r',pst2,sim_s,'b+');
lgd=legend('simulation data','theoretical results','predictive results');
lgd.FontSize = 20;
grid on;
xlabel('H','FontSize',20);
ylabel('MSE','FontSize',20);
title('\sigma(H)','FontSize',20);
saveas(gcf,sprintf('%s/res/pol_s_%d.png',root_path,sh_p));

