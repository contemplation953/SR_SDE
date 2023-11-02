function res=global_cost(coeff,pdf_matrix,H_list,D_list,moment1_matrix,moment2_matrix)
    %损失函数

    pdf_ans_matrix=zeros(size(pdf_matrix));
    monent1_ans_matrix=zeros(size(moment1_matrix));
    monent2_ans_matrix=zeros(size(moment1_matrix));
    c=coeff;

    for i=1:length(D_list)
        %d->x;h->y;
        x=D_list(i);

        m_x=@(y) c(1)+c(2).*x+c(3).*y+c(4).*x.*y+c(5).*y.^2;
        s_x=@(y) c(6)+c(7).*x+c(8).*y+c(9).*x.^2+c(10).*x.*y+c(11).*y.^2+c(12).*x.^2.*y+c(13).*x.*y.^2+c(14).*y.^3;
%         m_x=@(y) c(1)*x+c(2)*y.^2+c(3)*y+c(4);
%         s_x=@(y) c(5)*x.*y+c(6)*x.^2.*y.^2+c(7);
        [~, p_analytical]=fp_solve(H_list,m_x,s_x);
        
        pdf_ans_matrix(i,:)=p_analytical;
        monent1_ans_matrix(i,:)=m_x(H_list);
        stm=s_x(H_list);
        monent2_ans_matrix(i,:)=stm;
    end
    md=moment_displace(moment1_matrix,moment2_matrix,monent1_ans_matrix,monent2_ans_matrix);
    kld=KL_Divergence(pdf_matrix,pdf_ans_matrix,H_list);
    res=kld+md;
end

function res=KL_Divergence(p_data,p_analytical,pst)
    p_data(p_data<1e-20)=1e-20;
    p_analytical(p_analytical<1e-20)=1e-20;
    y=trapz(pst,p_data.*log(p_data./p_analytical),2);
    res=sum(y);
end

function res=moment_displace(moment1_matrix,moment2_matrix,monent1_ans_matrix,monent2_ans_matrix)
    d1=(moment1_matrix-monent1_ans_matrix).*(moment1_matrix-monent1_ans_matrix);
    d2=(moment2_matrix-monent2_ans_matrix).*(moment2_matrix-monent2_ans_matrix);
    res=sqrt(sum(sum(d1+d2)))/size(moment1_matrix,2);
end
