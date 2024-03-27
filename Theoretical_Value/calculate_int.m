function [M_int_list, Sigma_square_int_list, TT_int_list] = calculate_int(diff_list, R_table, data_path)
    len = length(diff_list);
    M_int_list = zeros(1,len);
    Sigma_square_int_list = zeros(1,len);
    TT_int_list = zeros(1,len);
    
    data_file=sprintf('%s/int_list.mat',data_path);
    if ~exist(data_file,'file')
	    data = load(data_file);
        M_int_list = data.M_int_list;
        Sigma_square_int_list = data.Sigma_square_int_list;
        TT_int_list = data.TT_int_list;
    else
        for H_index = 1: length(diff_list)
            H = diff_list(H_index);
            TT=integrate_gaussian(@TT_int, H, H_index, R_table);
            M_H=integrate_gaussian(@M_int, H, H_index, R_table);
            Sigma_H_square=integrate_gaussian(@Sigma_int, H, H_index, R_table);
            M_int_list(H_index) = M_H;
            Sigma_square_int_list(H_index) = Sigma_H_square;
            TT_int_list(H_index) = TT;
        end
        save(data_file,'M_int_list','Sigma_square_int_list', 'TT_int_list');
    end
end


function res = M_int(v1,v2,v3,v4,r,H)
    F=fun_F(v1,v2,v3,v4,r,H);
    eta1=fun_Eta(v1,v2,v3,v4);
    res=F.^(5/2).*(r.^2.*eta1+0.1).*r.^4.*sin(v1).^3.*sin(v2).^2.*sin(v3);
end

function res = Sigma_int(v1,v2,v3,v4,r,H)
    F=fun_F(v1,v2,v3,v4,r,H);
    res=F.^(5/2).*r.^4.*sin(v1).^3.*sin(v2).^2.*sin(v3);
end

function res = TT_int(v1,v2,v3,v4,r,H)
    F=fun_F(v1,v2,v3,v4,r,H);
    res=F.^(3/2).*r.^4.*sin(v1).^3.*sin(v2).^2.*sin(v3);
end


%%定义一些基础的函数
function res=fun_F(v1,v2,v3,v4,r,H)
    res=2.*H-fun_PhiB(v1,v2,v3,v4).*r.^2-fun_PhiA(v1,v2,v3,v4).*r.^4;
end

function res=fun_PhiA(v1,v2,v3,v4)
    res=(cos(v1)-sin(v1).*cos(v2)).^4/2 ...
    +(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).^4/2 ...
    +(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).^4/2 ...
    +(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).^4/2;
end

function res=fun_PhiB(v1,v2,v3,v4)
    res=1+(cos(v1)-sin(v1).*cos(v2)).^2 ...
    +(sin(v1).*cos(v2)-sin(v1).*sin(v2).*cos(v3)).^2 ...
    +(sin(v1).*sin(v2).*(cos(v3)-sin(v3).*cos(v4))).^2 ...
    +(sin(v1).*sin(v2).*sin(v3).*(cos(v4)-sin(v4))).^2;
end

function res=fun_Eta(v1,v2,v3,v4)
    res=0.01.*cos(v1).^2+0.01.*sin(v1).^2.*cos(v2).^2+0.01.*sin(v1).^2.*sin(v2).^2.*cos(v3).^2 ...
    +0.01.*sin(v1).^2.*sin(v2).^2.*sin(v3).^2.*cos(v4).^2+0.01.*sin(v1).^2.*sin(v2).^2.*sin(v3).^2.*sin(v4).^2;
end