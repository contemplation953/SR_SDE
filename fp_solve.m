function [pst,pdf] = fp_solve(pst, m, s)
    %solve 1d sde's FPK function steady pdf
    %dx=m(x)dt+s(x)dB(t)
    %p(x)=C/s^2(x)*exp(int(2*m(x)/s^2(x)))

    %%m(x) 漂移
    %%s(x) 扩散
    %pst bin


    fun=@(x) 2*m(x)./(s(x).^2);
    sigma_2x=@(x) (s(x).^2);

    pdf_list=zeros(1,length(pst));
    for i=1:length(pst)
        y=ComSimpson(fun,10,0,pst(i));
        pdf_list(i)=exp(y)/sigma_2x(pst(i));
    end
    pdf_list(pdf_list<0)=0;
    C=trapz(pst,pdf_list);

    pdf=pdf_list/C;
end