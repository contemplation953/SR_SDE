clc;clear;

CoreNum=25;
if isempty(gcp('nocreate'))
   parpool(CoreNum);
end

pm_matrix=zeros(2,25);
index=1;
for mi=1:5
    for si=1:5
        pm_matrix(1,index)=mi;
        pm_matrix(2,index)=si;
        index=index+1;
    end
end

root_path=pwd();

parfor i=1:25
    mi=pm_matrix(1,i);
    si=pm_matrix(2,i);
    optimize(mi,si,root_path);
end

delete(gcp('nocreate'));

function res=optimize(pm,ps,root_path)
    tic;
    step=1;
    val_list=zeros(length(1:step:61),1);
    x_list=zeros(61,pm+ps+2);
    for i=1:step:61
        data=load(sprintf('%s/data/H_%d.mat',root_path,i));
        h=0.01;
        H_list1=0.01:h:25-h;
        
        moment1=data.moment1(1:length(H_list1));
        moment2=sqrt(data.moment2(1:length(H_list1))); 
        coeff=init(H_list1,moment1,moment2,pm,ps);

        H_list=0.01:h:20-h;
        y1=polyval(coeff(1:pm+1),H_list);
        y2=polyval(coeff(pm+2:end),H_list);    
        p_data=data.p_data(1:length(H_list));
       
        options = optimset('MaxFunEvals',1e9,'MaxIter',1e3,'TolFun',1e-4);
        fun=@(coeff) ssr_cost(coeff,pm,ps,p_data,H_list,y1,y2);

        [x,val]=fminsearch(fun,coeff,options);
        x_list(i,:)=x;
        val_list(i)=val; 
    end
    data_path=sprintf('%s/data/SSR_%d_%d.mat',root_path,pm,ps);
    res=mean(val_list);
    save(data_path,'val_list','x_list','res','-v7.3');
    toc;
end


function coeff=init(pst,m1_data,m2_data,pm,ps)
    coeff1=polyfit(pst,m1_data,pm);
    coeff2=polyfit(pst,m2_data,ps);
    coeff=[coeff1,coeff2];
end

