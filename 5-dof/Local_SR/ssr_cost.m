function res=ssr_cost(coeff,pm,ps,p_data,pst,m1_data,m2_data)
    coeff1=coeff(1:pm+1);
    coeff2=coeff(pm+2:end);
    [m_x,s_x]=funcFactory(coeff1,coeff2);

    md=moment_displace(pst,m_x,s_x,m1_data,m2_data);
    [~, p_analytical]=fp_solve(pst,m_x,s_x);
    kld=KL_Divergence(p_data,p_analytical,pst);
    res=md+kld;
end

function res=KL_Divergence(p_data,p_analytical,pst)
    p_data(p_data<1e-20)=1e-20;
    p_analytical(p_analytical<1e-20)=1e-20;
    
    res=trapz(pst,p_data.*log(p_data./p_analytical));
end

function res=moment_displace(pst,m_x,s_x,m1_data,m2_data)
    m1_ans=m_x(pst);
    m2_ans=s_x(pst);
    d1=(m1_ans-m1_data);
    d2=(m2_ans-m2_data);
    
    res=sqrt(sum(d1.*d1+d2.*d2))/length(m1_data);
end


function [mh,sh]=funcFactory(coeff1,coeff2)
    p1=length(coeff1);
    p2=length(coeff2);
    
    str1=string(zeros(1,p1+1));
    str2=string(zeros(1,p2+1));

    str1(1)='@(x) ';
    str2(1)='@(x) ';
    for i=1:p1
        if i==p1
            str1(i+1)=[num2str(coeff1(i))];
        else
            str1(i+1)=[num2str(coeff1(i)),sprintf('*x.^%d+',p1-i)];
        end
    end

    for i=1:p2
        if i==p2
            str2(i+1)=[num2str(coeff2(i))];
        else
            str2(i+1)=[num2str(coeff2(i)),sprintf('*x.^%d+',p2-i)];
        end
    end
    str1=join(str1);
    str2=join(str2);

    mh=matlabFunction(str2sym(str1));
    sh=matlabFunction(str2sym(str2));
end