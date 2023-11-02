function R_table = calculate_R(n, H_list, data_path)
    
    data_file=sprintf('%s/R_table.mat',data_path);

    if exist(data_file,'file')
	    data = load(data_file);
        R_table = data.R_table;
    else
        root_point = getGQ_table(n, 'Legendre');
        R_table = calculate(n, root_point, H_list);
        save(data_file,'R_table');
    end
end

function res = calculate(n, root_point, H_list)
    syms A B R H;
    eq = A*R^4+B*R^2-2*H;
    symbol_res=solve(eq,R);
    symbol_fun=matlabFunction(symbol_res,'Vars',[A B H]);

    res = zeros(n, n, n, n, length(H_list));
    v4_inf=0;
    v4_sup=2*pi;
    v3_inf=0;
    v3_sup=pi;
    v2_inf=0;
    v2_sup=pi;
    v1_inf=0;
    v1_sup=pi;

    h1 = (v4_sup-v4_inf)/2;
    h2 = (v4_sup+v4_inf)/2;
    for i4 = 1 : n
        v4 = h1*root_point(i4)+h2;
        v4_1=(v3_sup-v3_inf)/2;
        v4_2=(v3_sup+v3_inf)/2;
    
        for i3 = 1 : n
            v3 = v4_1*root_point(i3)+v4_2;
            v3_1=(v2_sup-v2_inf)/2;
            v3_2=(v2_sup+v2_inf)/2;
    
            for i2 = 1 : n
                v2 = v3_1*root_point(i2)+v3_2;
                v2_1=(v1_sup-v1_inf)/2;
                v2_2=(v1_sup+v1_inf)/2;
                
                tic
                for i1 = 1 : n
                    v1 = v2_1*root_point(i2)+v2_2;
                    onevec=ones(1,length(H_list));
                    A_p=getAparam(v1*onevec,v2*onevec,v3*onevec,v4*onevec);
                    B_p=getBparam(v1*onevec,v2*onevec,v3*onevec,v4*onevec);
                    fun_val=max_real(symbol_fun(A_p,B_p,H_list));
                    for i = 1:length(fun_val)
                        res(i4, i3, i2, i1, i, 1) = fun_val(i);
                    end
                end
            end
        end
    end
end

function res = getAparam(v1,v2,v3,v4)
    %大概需要0.02s
    res=(cos(v1)-sin(v1).*cos(v2)).*(cos(v1)-sin(v1).*cos(v2)).*(cos(v1)-sin(v1).*cos(v2)).*(cos(v1)-sin(v1).*cos(v2))/2 ...
    +(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).*(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).*(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).*(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3))/2 ...
    +(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).*(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).*(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).*(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4)))/2 ...
    +(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).*(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).*(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).*(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4)))/2;


end

function res = getBparam(v1,v2,v3,v4)
    res=1+(cos(v1)-sin(v1).*cos(v2)).*(cos(v1)-sin(v1).*cos(v2)) ...
    +(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).*(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)) ...
    +(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).*(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))) ...
    +(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).*(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4)));
end

function res=max_real(data)
    res=zeros(1,size(data,2));
    for i=1:size(data,2)
        col=data(:,i);
        res(i)=max(col(col==real(col)));
    end
end
