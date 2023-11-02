clc, clear;

%10dof m2 s5
%20dof m2 s3
mh_p=2;
sh_p=5;
root_path=pwd();

data_index=1;
data1=load(sprintf('%s/N-DOF/10_dof_data/H_%d.mat',root_path,data_index));
pdf_data=data1.p_data;

data=load(sprintf('%s/N-DOF/10_dof_data/SSR_%d_%d.mat',root_path,mh_p,sh_p));

coeff=data.x_list(data_index,:);
coeff1=coeff(1:mh_p+1);
coeff2=coeff(mh_p+2:end);
[drift,diffusion]=FuncFactory(coeff1,coeff2);


%统计间隔
pst=0.01:0.01:30;
statistics_step=pst(2)-pst(1);
Stepsize=5e-3;
%运行总时间
FinalTime=3500;

%开始存储的时间
BeginTime = 2000;
%样本数
NumberOfSample=200;
NumberOfSubinterval = ceil(FinalTime/Stepsize);
begin_index = BeginTime/Stepsize;


statistics_count=zeros(length(pst),1);
InitialValue=1;
H_data=InitialValue*rand(NumberOfSample,1);


for i=1:(NumberOfSubinterval-1)
    dB = sqrt(Stepsize)*randn(NumberOfSample, 1);
    H_data=H_data+Stepsize*drift(H_data)+dB.*diffusion(H_data);
    if i>begin_index
        statistics_count = statistics_count+count(pst,H_data);
    end
end

%计算pdf
total_num=NumberOfSample*(FinalTime-BeginTime)/Stepsize;
pdf_ans=statistics_count/total_num/statistics_step;
pdf_ans=pdf_ans';
kl=KL_Divergence(pdf_data,pdf_ans,pst);
plot(pst,pdf_data,'b',pst,pdf_ans,'r');

title(sprintf('kl=%.4f',kl));


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

function res=KL_Divergence(p_data,p_analytical,pst)
    p_data(p_data==0)=1e-100;
    p_analytical(p_analytical==0)=1e-100;
    res=trapz(pst,p_data.*log(p_data./p_analytical),2);
end