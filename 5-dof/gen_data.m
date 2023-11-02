clc, clear;

data_path=sprintf('%s/5-dof/data',pwd);

param_struct.beta1=0.01;
param_struct.beta2=0.01;
param_struct.beta3=0.01;
param_struct.beta4=0.01;
param_struct.beta5=0.01;

param_struct.c1=0.02;
param_struct.c2=0.02;
param_struct.c3=0.02;
param_struct.c4=0.02;
param_struct.c5=0.02;

param_struct.omega1=1;
param_struct.omega2=1;
param_struct.omega3=1;
param_struct.omega4=1;
param_struct.omega5=1;

param_struct.a12=1;
param_struct.a23=1;
param_struct.a34=1;
param_struct.a45=1;

param_struct.b12=1;
param_struct.b23=1;
param_struct.b34=1;
param_struct.b45=1;

%-----------------------------------------------
omega1=param_struct.omega1;
omega2=param_struct.omega2;
omega3=param_struct.omega3;
omega4=param_struct.omega4;
omega5=param_struct.omega5;

a12=param_struct.a12;
a23=param_struct.a23;
a34=param_struct.a34;
a45=param_struct.a45;


b12=param_struct.b12;
b23=param_struct.b23;
b34=param_struct.b34;
b45=param_struct.b45;

%ͳ�Ƽ��
pst=0.01:0.01:25;
statistics_step=pst(2)-pst(1);
statistics_bound=pst(end)/statistics_step;
Stepsize=5e-3;
%������ʱ��
FinalTime=3500;

%��ʼ�洢��ʱ��
BeginTime=2000;
%������
NumberOfSample=200;
NumberOfSubinterval = ceil(FinalTime/Stepsize);
begin_index=BeginTime/Stepsize;

file_index=1;

for D=0.005:0.0005:0.035
tic;
    %ֱ������
    statistics_count=zeros(length(pst),1);
    %ͳ��ÿ��������ֱ���ڵ�����
    monment_count=zeros(length(pst),1);
    %һ�׵����ؼ���
    moment1_diff=zeros(length(pst),1);
    %���׵����ؼ���
    moment2_diff=zeros(length(pst),1);

    
    InitialValue=3;
    %�˿̵�ֵ
    q1 = InitialValue*rand(NumberOfSample,1);
    q2 = InitialValue*rand(NumberOfSample,1); 
    q3 = InitialValue*rand(NumberOfSample,1);
    q4 = InitialValue*rand(NumberOfSample,1); 
    q5 = InitialValue*rand(NumberOfSample,1);
    
    p1 = InitialValue*rand(NumberOfSample,1); 
    p2 = InitialValue*rand(NumberOfSample,1);
    p3 = InitialValue*rand(NumberOfSample,1); 
    p4 = InitialValue*rand(NumberOfSample,1);
    p5 = InitialValue*rand(NumberOfSample,1); 
    
    for i=1:(NumberOfSubinterval-1)
        %%-------------------------------------------------------------------
        %%�������εĴ���ÿһ�е�����Ҫ��������
        %%-------------------------------------------------------------------
        GW = sqrt(2*D/Stepsize)*randn(NumberOfSample, 5);
        yn = [q1,q2,q3,q4,q5,p1,p2,p3,p4,p5];
    
        k1 = Diff_Drift(yn,GW,param_struct);
        k2 = Diff_Drift(yn+Stepsize/2*k1,GW,param_struct);
        k3 = Diff_Drift(yn+Stepsize/2*k2,GW,param_struct);
        k4 = Diff_Drift(yn+Stepsize*k3,GW,param_struct);
        y_next = yn+(k1+2*k2+2*k3+k4)/6*Stepsize;

         %���㵱ǰ��������
        H = 1/2.*(p1.*p1+p2.*p2+p3.*p3+p4.*p4+p5.*p5) ...
            + 1/2.*(omega1.*omega1.*q1.*q1+omega2.*omega2.*q2.*q2+omega3.*omega3.*q3.*q3+omega4.*omega4.*q4.*q4+omega5.*omega5.*q5.*q5)...
            + 1/2.*(a12.*(q1-q2).*(q1-q2)+a23.*(q2-q3).*(q2-q3)+a34.*(q3-q4).*(q3-q4)+a45.*(q4-q5).*(q4-q5)) ...
            + 1/4.*(b12.*(q1-q2).*(q1-q2).*(q1-q2).*(q1-q2)+b23.*(q2-q3).*(q2-q3).*(q2-q3).*(q2-q3)+b34.*(q3-q4).*(q3-q4).*(q3-q4).*(q3-q4)+b45.*(q4-q5).*(q4-q5).*(q4-q5).*(q4-q5));
        monment_count=monment_count+count(pst,H);
        

        %������һ��������
        q1=y_next(:,1);q2=y_next(:,2);q3=y_next(:,3);q4=y_next(:,4);q5=y_next(:,5);
        p1=y_next(:,6);p2=y_next(:,7);p3=y_next(:,8);p4=y_next(:,9);p5=y_next(:,10);
       
        H_tau = 1/2.*(p1.*p1+p2.*p2+p3.*p3+p4.*p4+p5.*p5) ...
            + 1/2.*(omega1.*omega1.*q1.*q1+omega2.*omega2.*q2.*q2+omega3.*omega3.*q3.*q3+omega4.*omega4.*q4.*q4+omega5.*omega5.*q5.*q5)...
            + 1/2.*(a12.*(q1-q2).*(q1-q2)+a23.*(q2-q3).*(q2-q3)+a34.*(q3-q4).*(q3-q4)+a45.*(q4-q5).*(q4-q5)) ...
            + 1/4.*(b12.*(q1-q2).*(q1-q2).*(q1 -q2).*(q1-q2)+b23.*(q2-q3).*(q2-q3).*(q2-q3).*(q2-q3)+b34.*(q3-q4).*(q3-q4).*(q3-q4).*(q3-q4)+b45.*(q4-q5).*(q4-q5).*(q4-q5).*(q4-q5));
       
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
    
    %����pdf
    total_num=NumberOfSample*(FinalTime-BeginTime)/Stepsize;
    p_data=statistics_count/total_num/statistics_step;
    p_data=p_data';

    %���㵼����(�ų�0�ĸ���)
    monment_count(monment_count==0)=1;
    moment1=(moment1_diff./monment_count)';
    moment2=(moment2_diff./monment_count)';
    save(sprintf('%s/H_%d.mat',data_path,file_index),'p_data','moment1','moment2','monment_count','-v7.3');
    file_index=file_index+1;
toc;
end

function dy = Diff_Drift(yn,w,param_struct)
    dy = zeros(size(yn,1),size(yn,2));
    beta1=param_struct.beta1;
    beta2=param_struct.beta2;
    beta3=param_struct.beta3;
    beta4=param_struct.beta4;
    beta5=param_struct.beta5;
    
    c1=param_struct.c1;
    c2=param_struct.c2;
    c3=param_struct.c3;
    c4=param_struct.c4;
    c5=param_struct.c5;
    
    omega1=param_struct.omega1;
    omega2=param_struct.omega2;
    omega3=param_struct.omega3;
    omega4=param_struct.omega4;
    omega5=param_struct.omega5;
    
    a12=param_struct.a12;
    a23=param_struct.a23;
    a34=param_struct.a34;
    a45=param_struct.a45;
    
    
    b12=param_struct.b12;
    b23=param_struct.b23;
    b34=param_struct.b34;
    b45=param_struct.b45;

    q1 = yn(:,1);
    q2 = yn(:,2);
    q3 = yn(:,3);
    q4 = yn(:,4);
    q5 = yn(:,5);

    p1 = yn(:,6);
    p2 = yn(:,7);
    p3 = yn(:,8);
    p4 = yn(:,9);
    p5 = yn(:,10);

    
    dy(:,1)=p1;
    dy(:,2)=p2;
    dy(:,3)=p3;
    dy(:,4)=p4;
    dy(:,5)=p5;
    dy(:,6)=-(beta1.*q1.*q1+c1).*p1-(omega1*omega1.*q1+a12*(q1-q2)+b12.*(q1-q2).*(q1-q2).*(q1-q2))+w(:,1);
    dy(:,7)=-(beta2.*q2.*q2+c2).*p2-(omega2*omega2.*q2+a12*(q2-q1)+b12.*(q2-q1).*(q2-q1).*(q2-q1)+a23.*(q2-q3)+b23.*(q2-q3).*(q2-q3).*(q2-q3))+w(:,2);
    dy(:,8)=-(beta3.*q3.*q3+c3).*p3-(omega3*omega3.*q3+a23*(q3-q2)+b23.*(q3-q2).*(q3-q2).*(q3-q2)+a34.*(q3-q4)+b34.*(q3-q4).*(q3-q4).*(q3-q4))+w(:,3);
    dy(:,9)=-(beta4.*q4.*q4+c4).*p4-(omega4*omega4.*q4+a34*(q4-q3)+b34.*(q4-q3).*(q4-q3).*(q4-q3)+a45.*(q4-q5)+b45.*(q4-q5).*(q4-q5).*(q4-q5))+w(:,4);
    dy(:,10)=-(beta5.*q5.*q5+c5).*p5-(omega5*omega5.*q5+a45*(q5-q4)+b45.*(q5-q4).*(q5-q4).*(q5-q4))+w(:,5);
end

%ֱ���ֲ�
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