function res = ComSimpson(f,n,a,b)
    format long;
    h = (b-a)/n;
    d = f(a);
    for i = a+h:h:b-h 
        d = d + (2 * f(i));
    end
    for i = a+h/2:h:b-h/2 
        d = d + (4 * f(i));
    end
    d = d + f(b);
    res = (d * h / 6);