clc, clear;


%训练集
D_list=0.005:0.0005:0.035;
data_path=sprintf('%s/N-DOF/20_dof_res',pwd);

for data_index=1:1:length(D_list)
    tic;
    data=load(sprintf('%s/global_params.mat',data_path));
    c=data.x;
    d=D_list(data_index);
    
    drift=@(h) c(1)+c(2)*d+c(3)*h+c(4)*d.^2+c(5)*d.*h+c(6)*h.^2;
    diffusion=@(h) c(7)+c(8)*d+c(9)*h+c(10)*d.^2+c(11)*d.*h+c(12)*h.^2+c(13)*d.^2.*h+c(14)*d.*h.^2+c(15)*h.^3;

    
    %统计间隔
    pst=0.01:0.01:30;
    statistics_step=pst(2)-pst(1);
    Stepsize=5e-3;
    %运行总时间
    FinalTime=3500;
    
    %开始存储的时间
    BeginTime = 2000;
    %样本数
    NumberOfSample=800;
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

    save(sprintf('%s/H_%d.mat',data_path,data_index),'pdf_ans');
    
    toc;
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
