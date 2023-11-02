function num_res = integrate_gaussian(fun, H, H_index, R_table)

    if ~isa(fun,'function_handle') || nargin(fun) == 0
        error('integralN:BadIntegrand', ...
            ['First argument must be a function handle accepting one or more inputs.\n', ...
            'To integrate a constant c, integrate @(x,y)c*ones(size(x)) for 2D,\n', ...
            '@(x,y,z)c*ones(size(x)) for 3D, etc.']);
    end
    format long;
    glq_n = 5;
    [vloop4,vloop3,vloop2,vloop1,rloop]=deal(glq_n);
    [root_point,weight] = getGQ_table(glq_n, 'Legendre');
    
    
    v4_inf=0;
    v4_sup=2*pi;
    v3_inf=0;
    v3_sup=pi;
    v2_inf=0;
    v2_sup=pi;
    v1_inf=0;
    v1_sup=pi;
    r_inf=0;
    
    h1 = (v4_sup-v4_inf)/2;
    h2 = (v4_sup+v4_inf)/2;
    J=0;
    
    for i4 = 1 : vloop4
        JX_4 = 0;
        v4 = h1*root_point(i4)+h2;
        sup_v=v3_sup;
        inf_v=v3_inf;
        v4_1=(sup_v-inf_v)/2;
        v4_2=(sup_v+inf_v)/2;
    
        for i3 = 1 : vloop3
            JX_3 = 0;
            v3 = v4_1*root_point(i3)+v4_2;
            sup_v=v2_sup;
            inf_v=v2_inf;
            v3_1=(sup_v-inf_v)/2;
            v3_2=(sup_v+inf_v)/2;
    
            for i2 = 1 : vloop2
                JX_2 = 0;
                v2 = v3_1*root_point(i2)+v3_2;
                sup_v=v1_sup;
                inf_v=v1_inf;
                v2_1=(sup_v-inf_v)/2;
                v2_2=(sup_v+inf_v)/2;
    
                for i1 = 1 : vloop1
                    JX_1 = 0;
                    v1 = v2_1*root_point(i2)+v2_2;
                    sup_v=getR(i4,i3,i2,i1,H_index,R_table);
                    inf_v=r_inf;
                    v1_1=(sup_v-inf_v)/2;
                    v1_2=(sup_v+inf_v)/2;
    
                    for i = 1 : rloop
                        r_val = v1_1*root_point(i)+v1_2;
                        Q = fun(v1,v2,v3,v4,r_val,H);
                        JX_1 = JX_1 + weight(i)*Q;
                    end
                    JX_2 = JX_2 + weight(i1)*v1_1*JX_1;
                end
                JX_3 = JX_3 + weight(i2)*v2_1*JX_2;
            end
            JX_4 = JX_4 + weight(i3)*v3_1*JX_3;
        end
        J = J + weight(i4)*v4_1*JX_4;
    end
    num_res = h1*J;
end

function r_val = getR(i4, i3, i2, i1, H_index, R_table)
    r_val = R_table(i4, i3, i2, i1, H_index, 1);
end
