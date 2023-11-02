clc;clear;


%训练集
DD_list=0.005:0.0005:0.035;
%检验集
VD_list=[0.005,0.015,0.025,0.035];
%统计间隔
H_list=0.01:0.01:25;

D_matrix=repmat(DD_list',1,length(H_list));
H_matrix=repmat(H_list,length(DD_list),1);
pdf_data_matrix=zeros(length(DD_list),length(H_list));

root_path=pwd();
data_path=sprintf('%s/5-dof/data',pwd);

for i=1:length(DD_list)
    tic;
    data1=load(sprintf('%s/H_%d.mat',data_path,i));
    pdf_data_matrix(i,:)=data1.p_data;
end

pdf_data_matrix=zeros(length(VD_list),length(H_list));
for i=1:length(VD_list)
    tic;
    data1=load(sprintf('%s/val_H_%d.mat',data_path,i));
    pdf_data_matrix(i,:)=data1.p_data;
end

figure(1);
plot(H_list,pdf_data_matrix(1,:),'r-','LineWidth',2);
hold on;
plot(H_list,pdf_data_matrix(2,:),'g-','LineWidth',2);
hold on;
plot(H_list,pdf_data_matrix(3,:),'b-','LineWidth',2);
hold on;
plot(H_list,pdf_data_matrix(4,:),'k-','LineWidth',2);

xlabel('H');
ylabel('P(H)');
legend('D=0.004','D=0.015','D=0.025','D=0.035');
hold off;

saveas(gcf,sprintf('%s/res-global/pdf_val.png',root_path));
