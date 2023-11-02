function [mh,sh]=FuncFactory(coeff1,coeff2)
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