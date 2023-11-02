clc;clear;

pm_matrix=zeros(2,25);
index=1;
for mi=2:5
    for si=1:5
        pm_matrix(1,index)=mi;
        pm_matrix(2,index)=si;
        index=index+1;
    end
end

data_path=sprintf('%s/N-DOF/20_dof_data',pwd());
val_list=zeros(25,1);

for i=1:25
    mi=pm_matrix(1,i);
    si=pm_matrix(2,i);
    optimize(mi,si,data_path);
end


function res=optimize(pm,ps,data_path)
    tic;
    step=10;
    val_list=zeros(length(1:step:61),1);
    x_list=zeros(61,pm+ps+2);
    for i=1:step:61
        data=load(sprintf('%s/H_%d.mat',data_path,i));
        H_list=0.01:0.01:30;
    
        p_data=data.p_data(1:length(H_list));
        moment1=data.moment1(1:length(H_list));
        moment2=data.moment2(1:length(H_list)); 

        pdf_bound=1e-3;
        idnex=find(p_data>=pdf_bound);
        begin_idnex=idnex(1);
        end_index=idnex(end);
        H_list_p=H_list(begin_idnex:end_index);
        moment1_p=moment1(begin_idnex:end_index);
        moment2_p=moment2(begin_idnex:end_index);

        coeff=init(H_list_p,moment1_p,moment2_p,pm,ps);
        options = optimset('MaxFunEvals',1e9,'MaxIter',1e3,'TolFun',1e-4);
        fun=@(coeff) nd_cost(coeff,pm,ps,p_data,H_list,moment1,moment2);

        [x,val]=fminsearch(fun,coeff,options);
        x_list(i,:)=x;
        val_list(i)=val; 
        val
    end
    save_path=sprintf('%s/SSR_%d_%d.mat',data_path,pm,ps);
    res=mean(val_list);
    save(save_path,'val_list','x_list','res','-v7.3');
    toc;
end


function coeff=init(pst,m1_data,m2_data,pm,ps)
    coeff1=polyfit(pst,m1_data,pm);
    coeff2=polyfit(pst,sqrt(m2_data),ps);
    coeff=[coeff1,coeff2];
end

