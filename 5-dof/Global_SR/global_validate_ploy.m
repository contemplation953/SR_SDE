clc;clear;

%产生检验集 train:1;validate:2
Validata_Flag=2;

data_path=sprintf('%s/5-dof/data',pwd);
ans_path=sprintf('%s/Theoretical_Value/data',pwd);

%训练集
DD_list=0.005:0.0005:0.035;
%检验集
VD_list=[0.004,0.0258,0.036,0.04];

if Validata_Flag==1
    D_list=DD_list;
else
    D_list=VD_list;
end


figure('visible','off');
for data_index=1:length(D_list)
    if Validata_Flag==1
        data=load(sprintf('%s/H_%d.mat',data_path,data_index));
        data1=load(sprintf('%s/analyse_%d.mat',ans_path,data_index));
    else
        data=load(sprintf('%s/val_H_%d.mat',data_path,data_index));
        data1=load(sprintf('%s/val_analyse_%d.mat',ans_path,data_index));
    end

    
    pst=0.01:0.01:20;
    pdf_bound=1e-3;
    p_data=data.p_data;
    idnex=find(p_data>pdf_bound);
    begin_idnex=idnex(1);
    end_index=idnex(end);
    pst=pst(begin_idnex:end_index);
    moment1=data.moment1(begin_idnex:end_index);
    moment2=data.moment2(begin_idnex:end_index);
    
    pst1=0.001:0.001:20;
    M_H=data1.M_H(1:length(pst1));
    Sigma_H_square=data1.Sigma_H_square(1:length(pst1));
    
    data2=load(sprintf('%s/global_params.mat',data_path));
    c=data2.x;
    pst2=0.5:0.5:20;
    
    x=D_list(data_index);
    m_x=@(y) c(1)+c(2).*x+c(3).*y+c(4).*x.*y+c(5).*y.^2;
    s_x=@(y) (c(6)+c(7).*x+c(8).*y+c(9).*x.^2+c(10).*x.*y+c(11).*y.^2+c(12).*x.^2.*y+c(13).*x.*y.^2+c(14).*y.^3).^2;
    
    sim_m=m_x(pst2);
    sim_s=s_x(pst2);

    if Validata_Flag==1
        png_path1=sprintf('%s/res-global/%d_moment1.png',pwd,data_index);
        png_path2=sprintf('%s/res-global/%d_moment2.png',pwd,data_index);
    else
        png_path1=sprintf('%s/res-global/val_%d_moment1.png',pwd,data_index);
        png_path2=sprintf('%s/res-global/val_%d_moment2.png',pwd,data_index);
    end

    plot(pst,moment1,'kO',pst1,M_H,'r',pst2,sim_m,'b+');
    grid on;
    xlabel('H','FontSize',20);
    ylabel('MSE','FontSize',20);
    lgd=legend('simulation data','theoretical results','predictive results');
    lgd.FontSize = 20;
    title('m(H)','FontSize',20);
    saveas(gcf,png_path1);

    plot(pst,moment2,'kO',pst1,Sigma_H_square,'r',pst2,sim_s,'b+');
    grid on;
    xlabel('H','FontSize',20);
    ylabel('MSE','FontSize',20);
    lgd=legend('simulation data','theoretical results','predictive results');
    lgd.FontSize = 20;
    title('\sigma^2(H)','FontSize',20);
    saveas(gcf,png_path2);

    
end

figure('visible','on');
