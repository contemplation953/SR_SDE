clc;clear;

pst=0.01:0.01:20;
mh_val=zeros(61,length(pst));
sh_val=zeros(61,length(pst));
d_list=0.005:0.0005:0.035;
d_val=repmat(d_list',1,length(pst));
h_val=repmat(pst,61,1);

data_path=sprintf('%s/Theoretical_Value/data',pwd);
for i=1:61
    data1=load(sprintf('%s/analyse_%d.mat',data_path,i));
    M_H=data1.M_H;
    Sigma_H_square=data1.Sigma_H_square;
    mh_val(i,:)=M_H(1:10:end*2/3);
    sh_val(i,:)=sqrt(Sigma_H_square(1:10:end*2/3));
end

fd_val=reshape(d_val',1,numel(d_val))';
fh_val=reshape(h_val',1,numel(h_val))';
fmh_val=reshape(mh_val',1,numel(mh_val))';
fsh_val=reshape(sh_val',1,numel(sh_val))';

figure(1);
surf(d_val,h_val,mh_val);
xlabel('D');
ylabel('H');
zlabel('M(H)');
shading interp;
sf1=fit([fd_val, fh_val],fmh_val,'poly12');
params_m=coeffvalues(sf1);
figure(2);
plot(sf1,[fd_val, fh_val],fmh_val);


figure(3);
surf(d_val,h_val,sh_val);
xlabel('D');
ylabel('H');
zlabel('\sigma^2(H)');
shading interp;
sf2=fit([fd_val, fh_val],fsh_val,'poly22');
figure(4);
plot(sf2,[fd_val, fh_val],fsh_val);
params_s=coeffvalues(sf2);
if params_s(1)<0
    params_s(1)=0;
end

save(sprintf('%s/ML_ans_params.mat',data_path),'fd_val','fh_val','fmh_val','fsh_val','-v7.3');