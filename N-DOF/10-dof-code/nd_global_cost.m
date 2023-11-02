function res=nd_global_cost(coeff,pdf_matrix,H_list,D_list,moment1_matrix,moment2_matrix)
    %损失函数
    pdf_ans_matrix=zeros(size(pdf_matrix));
    monent1_ans_matrix=zeros(size(moment1_matrix));
    monent2_ans_matrix=zeros(size(moment1_matrix));
    c=coeff;

    for i=1:length(D_list)
        d=D_list(i);

        m_x=@(h) c(1)+c(2)*d+c(3)*h+c(4)*d.^2+c(5)*d.*h+c(6)*h.^2;
        s_x=@(h) c(7)+c(8)*d+c(9)*h+c(10)*d.^2+c(11)*d.*h+c(12)*h.^2+c(13)*d.^2.*h+c(14)*d.*h.^2+c(15)*h.^3;

        [~, p_analytical]=fp_solve_inf(H_list,m_x,s_x);
        if(~isempty(find(p_analytical<0, 1)))
            res=1e8-rand*1e2;
            return;
        end
        pdf_ans_matrix(i,:)=p_analytical;
        monent1_ans_matrix(i,:)=m_x(H_list);
        stm=s_x(H_list);
        monent2_ans_matrix(i,:)=stm.*stm;
    end
    md=moment_displace(moment1_matrix,moment2_matrix,monent1_ans_matrix,monent2_ans_matrix,pdf_matrix);
    kld=KL_Divergence(pdf_matrix,pdf_ans_matrix,H_list);
    res=kld+md;
end

function res=KL_Divergence(p_data,p_analytical,pst)
    p_data(p_data<1e-100)=1e-100;
    p_analytical(p_analytical<1e-100)=1e-100;
    y=trapz(pst,p_analytical.*log(p_analytical./p_data),2);
    res=sum(y);
end

function res=moment_displace(moment1_matrix,moment2_matrix,monent1_ans_matrix,monent2_ans_matrix,pdf_matrix)
    d1=moment1_matrix-monent1_ans_matrix;
    d2=moment2_matrix-monent2_ans_matrix;
    
    d11=0;
    d22=0;
    for i=1:size(d1,1)
        d11=d11+sqrt(d1(i,:).*d1(i,:))*pdf_matrix(i,:)';
        d22=d22+sqrt(d2(i,:).*d2(i,:))*pdf_matrix(i,:)';
    end
    res=(d11+d22)/numel(moment1_matrix);
end
