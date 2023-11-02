function [pst,pdf] = fp_solve_2(pst, m, s)
    %solve 1d sde's FPK function steady pdf
    %dx=m(x)dt+s(x)dB(t)
    %p(x)=C/s^2(x)*exp(int(2*m(x)/s^2(x)))

    %%m(x) 漂移
    %%s(x) 扩散项的平方
    %pst bin


    fun=@(x) 2*m(x)./s(x);

    pdf_list=zeros(1,length(pst));
    for i=1:length(pst)
        y=ComSimpson(fun,10,0,pst(i));
        pdf_list(i)=exp(y)/s(pst(i));
    end
    pdf_list(pdf_list<1e-20)=1e-20;
    C=trapz(pst,pdf_list);

    pdf=pdf_list/C;
end