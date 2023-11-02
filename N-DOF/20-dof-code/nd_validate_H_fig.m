clc;clear;

%训练集
DD_list=0.005:0.0005:0.035;
%检验集
VD_list=[0.005,0.015,0.025,0.35];
%统计间隔
H_list=0.01:0.01:60;
data_path=sprintf('%s/N-DOF',pwd);

D_matrix=repmat(DD_list',1,length(H_list));
H_matrix=repmat(H_list,length(DD_list),1);
pdf_data_matrix=zeros(length(DD_list),length(H_list));
pdf_mlres_matrix=zeros(length(DD_list),length(H_list));

tic;
for i=1:length(DD_list)
    data1=load(sprintf('D:/ChenXiYuan/workspace/SMl/sde/Ham-5-DOF/N-DOF/20_dof_data/H_%d.mat',i));
    pdf_data_matrix(i,:)=data1.p_data;
    
    data2=load(sprintf('D:/ChenXiYuan/workspace/SMl/sde/Ham-5-DOF/N-DOF/20_dof_res/H_%d.mat',i));
    pdf_mlres_matrix(i,:)=data2.pdf_ans;
end

figure(1);
surfl(H_matrix,D_matrix,pdf_data_matrix);
shading interp;
xlabel('H');
ylabel('D');
zlabel('pdf');
saveas(gcf,sprintf('%s/20_dof_res/20dof_pdf_3d_data.png',data_path));

figure(2);
surfl(H_matrix,D_matrix,pdf_mlres_matrix);
shading interp;
xlabel('H');
ylabel('D');
zlabel('pdf');
saveas(gcf,sprintf('%s/20_dof_res/20dof_pdf_3d_ml.png',data_path));

pdf_data_matrix=zeros(length(VD_list),length(H_list));
pdf_mlres_matrix=zeros(length(VD_list),length(H_list));
for i=1:length(VD_list)
    tic;
    data1=load(sprintf('D:/ChenXiYuan/workspace/SMl/sde/Ham-5-DOF/N-DOF/20_dof_data/H_%d.mat',i));
    pdf_data_matrix(i,:)=data1.p_data;
    
    data2=load(sprintf('D:/ChenXiYuan/workspace/SMl/sde/Ham-5-DOF/N-DOF/20_dof_res/H_%d.mat',i));
    pdf_mlres_matrix(i,:)=data2.pdf_ans;
end

figure(3);
plot(H_list(1:5:end/5),pdf_data_matrix(1,1:5:end/5),'rO',H_list(1:end/2),pdf_mlres_matrix(1,1:end/2),'r-');
hold on;
plot(H_list(1:40:end*0.5),pdf_data_matrix(2,1:40:end*0.5),'gO',H_list(1:end/2),pdf_mlres_matrix(2,1:end/2),'g-');
hold on;
plot(H_list(1:40:end*0.5),pdf_data_matrix(3,1:40:end*0.5),'bO',H_list(1:end/2),pdf_mlres_matrix(3,1:end/2),'b-');
hold on;
plot(H_list(1:40:end*0.5),pdf_data_matrix(4,1:40:end*0.5),'kO',H_list(1:end/2),pdf_mlres_matrix(4,1:end/2),'k-');
xlabel('H','FontSize',20);
ylabel('P(H)','FontSize',20);
legend('D=0.005 simulation data','D=0.005 predictive results','D=0.015 simulation data','D=0.015 predictive results','D=0.025 simulation data','D=0.025 predictive results','D=0.035 simulation data','D=0.035 predictive results');
hold off;
saveas(gcf,sprintf('%s/20_dof_res/20dof_pdf_demo.png',data_path));
