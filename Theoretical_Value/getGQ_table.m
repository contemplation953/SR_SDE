%获取Gaussian Legendre Quadrature 参数表
function [X,C] = getGQ_table(n, type)
    switch type
       case 'Legendre'
          [X,C] = getGauss_Legendre_table(n);
       case 'Laguerre'
          [X,C] = getGauss_Laguerre_table(n);
       otherwise
          X = 0;
          C = 0;
    end
end

function [X,C] = getGauss_Legendre_table(n)
    syms x;
    ft = ( x * x - 1)^n;
    y = diff( ft , n);
    str = prod(1 : n);
    f = (1 / (2^n * str))*y;
    f = collect(f);
    
    if n < 6
        X = eval(solve(f));
    else
        X = solve_polynomial(matlabFunction(f),[-1,1],n);
    end
    X = sort(X,'descend');
    C = zeros(1, length(X));
    for i = 1 : n
        int_val = 1;
        for j = 1 : n
            if i == j
                continue;
            else
                int_val = int_val*(x-X(j))/(X(i)-X(j));
            end
        end
        C(i) = integral(matlabFunction(int_val),-1,1);
    end
end


function [X,C] = getGauss_Laguerre_table(n)
    syms x;
    L_n = exp(x)*diff(exp(-x)*x^n,n);
    X = eval(solve(L_n));
    X = sort(X,'descend');

    dL_n=diff(L_n,1);
    
    C = zeros(1, length(X));
    fh = matlabFunction(dL_n);
    for i = 1 : n
        C(i) = prod(1 : n)^2/(X(i)*(fh(X(i)))^2);
    end
end
