clc;clear;


%训练集
DD_list=0.005:0.0005:0.035;
%检验集
VD_list=[0.004,0.0258,0.036];
%统计间隔
H_list=0.01:0.01:25;

D_matrix=repmat(DD_list',1,length(H_list));
H_matrix=repmat(H_list,length(DD_list),1);
pdf_data_matrix=zeros(length(DD_list),length(H_list));
pdf_mlres_matrix=zeros(length(DD_list),length(H_list));

root_path=pwd();
data_path=sprintf('%s/5-dof/data',pwd);

% for i=1:length(DD_list)
%     tic;
%     data1=load(sprintf('%s/H_%d.mat',data_path,i));
%     pdf_data_matrix(i,:)=data1.p_data;
%     
%     data2=load(sprintf('%s/res-global/data/H_%d.mat',root_path,i));
%     pdf_mlres_matrix(i,:)=data2.pdf_ans;
% end
% 
% figure(1);
% surfl(H_matrix,D_matrix,pdf_data_matrix);
% shading interp;
% xlabel('H');
% ylabel('D');
% zlabel('pdf');
% saveas(gcf,sprintf('%s/res-global/data_pdf_3d.png',root_path));
% 
% figure(2);
% surfl(H_matrix,D_matrix,pdf_mlres_matrix);
% shading interp;
% xlabel('H');
% ylabel('D');
% zlabel('pdf');
% saveas(gcf,sprintf('%s/res-global/mlres_pdf_3d.png',root_path));

pdf_data_matrix=zeros(length(VD_list),length(H_list));
pdf_mlres_matrix=zeros(length(VD_list),length(H_list));
for i=1:length(VD_list)
    tic;
    data1=load(sprintf('%s/val_H_%d.mat',data_path,i));
    pdf_data_matrix(i,:)=data1.p_data;
    
    data2=load(sprintf('%s/res-global/data/val_H_%d.mat',root_path,i));
    pdf_mlres_matrix(i,:)=data2.pdf_ans;
end

figure(3);
plot(H_list(1:3:end/10),pdf_data_matrix(1,1:3:end/10),'rO',H_list(1:end/8),pdf_mlres_matrix(1,1:end/8),'r-');
hold on;
plot(H_list(1:40:end*0.6),pdf_data_matrix(2,1:40:end*0.6),'gO',H_list(1:end*0.6),pdf_mlres_matrix(2,1:end*0.6),'g-');
hold on;
plot(H_list(1:40:end*0.8),pdf_data_matrix(3,1:40:end*0.8),'bO',H_list(1:end*0.8),pdf_mlres_matrix(3,1:end*0.8),'b-');
% hold on;
% plot(H_list(1:40:end),pdf_data_matrix(4,1:40:end),'kO',H_list,pdf_mlres_matrix(4,:),'k-');

xlabel('H');
ylabel('P(H)');
legend('D=0.004 simulation data','D=0.004 predictive results','D=0.0258 simulation data','D=0.0258 predictive results','D=0.036 simulation data','D=0.036 predictive results','D=0.04 simulation data','D=0.04 predictive results');
hold off;

% saveas(gcf,sprintf('%s/res-global/pdf_diff.png',root_path));
