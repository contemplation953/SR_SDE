clc, clear;

%N-dof system
dof_N=20;

%统计间隔
pst=0.01:0.01:60;
statistics_step=pst(2)-pst(1);
statistics_bound=pst(end)/statistics_step;
Stepsize=5e-3;
%运行总时间
FinalTime=350;

%开始存储的时间
BeginTime=200;
%样本数
NumberOfSample=200;
NumberOfSubinterval = ceil(FinalTime/Stepsize);
begin_index=BeginTime/Stepsize;

root_path=pwd();

D_list=0.005:0.0005:0.035;

for par_index=1:1
tic;
    D=D_list(par_index);
    %直方计数
    statistics_count=zeros(length(pst),1);
    %统计每个导数矩直方内的数量
    monment_count=zeros(length(pst),1);
    %一阶导数矩计数
    moment1_diff=zeros(length(pst),1);
    %二阶导数矩计数
    moment2_diff=zeros(length(pst),1);

    
    InitialValue=1;
    %此刻的值 1:N are q_i;N+1:end are p_i [q1,q2,...qN,p1,p2,...,pN]
    yn=InitialValue*rand(NumberOfSample,2*dof_N);
    
    for i=1:(NumberOfSubinterval-1)
        %计算当前步的能量
        H = cal_H(dof_N,yn);
        monment_count=monment_count+count(pst,H);

        %%-------------------------------------------------------------------
        %%极其隐蔽的错误，每一列的噪声要单独生成
        %%-------------------------------------------------------------------
        GW = sqrt(2*D/Stepsize)*randn(NumberOfSample, dof_N);
    
        k1 = diff_fun(dof_N,yn,GW);
        k2 = diff_fun(dof_N,yn+Stepsize/2*k1,GW);
        k3 = diff_fun(dof_N,yn+Stepsize/2*k2,GW);
        k4 = diff_fun(dof_N,yn+Stepsize*k3,GW);
        yn = yn+(k1+2*k2+2*k3+k4)/6*Stepsize;  

        %计算下一步的能量
        H_tau = cal_H(dof_N,yn);
       
        m1=(H_tau-H)/Stepsize;
        m2=(H_tau-H).*(H_tau-H)/Stepsize;
        h_index=ceil(H/statistics_step);
        h_index(h_index>statistics_bound)=statistics_bound;
        for l1=1:length(H)
            moment1_diff(h_index(l1))=moment1_diff(h_index(l1))+m1(l1);
            moment2_diff(h_index(l1))=moment2_diff(h_index(l1))+m2(l1);
        end
    
        if i > begin_index
            statistics_count = statistics_count+count(pst,H);
        end
        
    end
    
    %计算pdf
    total_num=NumberOfSample*(FinalTime-BeginTime)/Stepsize;
    p_data=statistics_count/total_num/statistics_step;
    p_data=p_data';

    %计算导数矩(排除0的干扰)
    monment_count(monment_count==0)=1;
    moment1=(moment1_diff./monment_count)';
    moment2=(moment2_diff./monment_count)';
    parsave(sprintf('%s/20_dof_data/H_%d.mat',root_path,par_index),p_data,moment1,moment2,monment_count);
toc;
end
plot(p_data);

function dy=diff_fun(N,yn,w)
    dy=zeros(size(yn,1),size(yn,2));
    q=yn(:,1:N);
    p=yn(:,N+1:end);

    for i=1:N
        fi=0.01*q(:,i).*q(:,i)+0.02;
        if i==1
            gi=q(:,1)+q(:,1)-q(:,2)+(q(:,1)-q(:,2)).*(q(:,1)-q(:,2)).*(q(:,1)-q(:,2));
        elseif i==N
            gi=q(:,N)+q(:,N)-q(:,N-1)+(q(:,N)-q(:,N-1)).*(q(:,N)-q(:,N-1)).*(q(:,N)-q(:,N-1));
        else
            gi=q(:,i)+q(:,i)-q(:,i-1)+(q(:,i)-q(:,i-1)).*(q(:,i)-q(:,i-1)).*(q(:,i)-q(:,i-1)) ...
                +q(:,i)-q(:,i+1)+(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1));
        end
        dy(:,i)=p(:,i);
        dy(:,N+i)=-fi.*p(:,i)-gi+w(:,i);
    end
end

function H=cal_H(N,yn)
    q=yn(:,1:N);
    p=yn(:,N+1:end);
    H=0;
    for i=1:N
        H=H+p(:,i).*p(:,i)+q(:,i).*q(:,i);
    end
    H=H/2;

    for i=1:N-1
        H=H+1/2*(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1)) ...
            +1/4*(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1)).*(q(:,i)-q(:,i+1));
    end
end


%直方分布
function statistics_count=count(pst,data)
    N=length(pst);
    statistics_count=zeros(N,1);
    range_data=data(data>pst(1) & data<pst(end));
    h=pst(2)-pst(1);
    index=ceil(range_data/h);
    for i=1:length(range_data)
        statistics_count(index(i))=statistics_count(index(i))+1;
    end
end

function parsave(fname, p_data,moment1,moment2,monment_count)
  save(fname,'p_data','moment1','moment2','monment_count','-v7.3');
end