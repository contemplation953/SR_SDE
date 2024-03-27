clear;clc;


data_path=sprintf('%s/Theoretical_Value/data',pwd);

%用于估算概率密度
step=0.1;
H_list = step:step:20;

%用于算积分
diff_step=0.001;
diff_list=diff_step:diff_step:30;

%提前计算出每个节点R的取值
R_table=calculate_R(5, diff_list, data_path);

[M_int_list, Sigma_square_int_list, TT_int_list] = calculate_int(diff_list, R_table, data_path);

file_index=1;

% D_add=[0.004,0.0258,0.036,0.04];
D_list=0.005:0.0005:0.035;
for D=D_list
    eta3=D*5;
    M_H=eta3 - 1./(5.*TT_int_list).*M_int_list;
    Sigma_H_square=2*eta3./(5.*TT_int_list).*Sigma_square_int_list;
    
    pdf_list = zeros(1,length(H_list));
    for i = 1 : length(H_list) 
        pdf_list(i) = getPDF(i * 100, M_int_list, Sigma_square_int_list, TT_int_list, D, diff_step);
    end
    C=trapz(H_list,pdf_list);
    pdf=pdf_list/C;
    %save(sprintf('%s/val_analyse_%d.mat',data_path,file_index),'M_H','Sigma_H_square','pdf','-v7.3');
    save(sprintf('%s/analyse_%d.mat',data_path,file_index),'M_H','Sigma_H_square','pdf','-v7.3');
    file_index=file_index+1;
end


function p_H = getPDF(index, M_int_list, Sigma_square_int_list, TT_int_list, D, diff_step)
    eta3=D*5;
    tt_val = TT_int_list(index);
    sigma2_val=2*eta3/(5*tt_val)*Sigma_square_int_list(index);


    TT = TT_int_list(1:index);
    M_H=eta3 - 1./(5.*TT(1:index)).*M_int_list(1:index);
    Sigma_H_square=2*eta3./(5.*TT(1:index)).*Sigma_square_int_list(1:index);
    
    int_f = 2*M_H./Sigma_H_square;
    int_sum = sum((int_f(1:index-1)+int_f(2:index))/2)*diff_step;
    
    p_H = exp(int_sum)/sigma2_val;
end
