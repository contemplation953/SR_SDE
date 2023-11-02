clc, clear;

root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);

data_index=1;
data1=load(sprintf('%s/H_%d.mat',data_path,data_index));
pdf_data=data1.p_data;

mh_p=4;
sh_p=5;
root_path=pwd();
data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));

coeff=data.x_list(data_index,:);
coeff1=coeff(1:mh_p+1);
coeff2=coeff(mh_p+2:end);
[drift,diffusion]=FuncFactory(coeff1,coeff2);


%统计间隔
pst=0.01:0.01:25;
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
figure(1);
plot(pst(1:2:end),pdf_data(1:2:end),'rO',pst,pdf_ans,'b');
xlabel('H','FontSize',20);
ylabel('PDF','FontSize',20);
lgd=legend('simulation data','predictive results');
lgd.FontSize = 20;
title(sprintf('D_{KL}=%.4f',kl),'FontSize',20);
saveas(gcf,sprintf('%s/pdf_%d_%d.png',ssr_path,mh_p,sh_p));

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
    p_data(p_data<1e-20)=1e-20;
    p_analytical(p_analytical<1e-20)=1e-20;
    res=trapz(pst,p_data.*log(p_data./p_analytical));
end