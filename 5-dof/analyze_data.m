clc;clear;

mh_p=4;
sh_p=5;
root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);
data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));

h=0.01;
H_list=0.01:h:20-h;
d_list=0.005:0.0005:0.035;

fd_max=zeros(61,length(H_list));
fh_max=zeros(61,length(H_list));
fmh_max=zeros(61,length(H_list));
fsh_max=zeros(61,length(H_list));
for i=1:61
    coeff=data.x_list(i,:);
    coeff1=coeff(1:mh_p+1);
    coeff2=coeff(mh_p+2:end);
    [drift,diffusion]=FuncFactory(coeff1,coeff2);

    sim_m=drift(H_list);
    sim_s=(diffusion(H_list));
    fmh_max(i,:)=sim_m';
    fsh_max(i,:)=sim_s';
    fd_max(i,:)=ones(1,length(H_list))*d_list(i);
    fh_max(i,:)=H_list;
end

fd_val=reshape(fd_max,numel(fd_max),1);
fh_val=reshape(fh_max,numel(fh_max),1);
fmh_val=reshape(fmh_max,numel(fmh_max),1);
fsh_val=reshape(fsh_max,numel(fsh_max),1);
%% 
figure(1);
surf(fd_max,fh_max,fmh_max);
shading interp;
grid on;
xlabel('D');
ylabel('H');
zlabel('M(H)');


% coff_init1=[-0.04005270118,0.04147547001,5.959692263518988,1,-5.959692263518988];
% fun_mh=@(c1,c2,c3,c4,c5,x,y) c1*y+c2+c3*x+c4*x.^2.*y+c5*x.^2;
% sf1=fit([fd_val, fh_val],fmh_val,fun_mh,'StartPoint',coff_init1);
sf1=fit([fd_val, fh_val],fmh_val,'poly12');
mhv1=reshape(sf1([fd_val, fh_val]),61,length(H_list));

params_m=coeffvalues(sf1);
figure(2);
plot(sf1,[fd_val, fh_val],fmh_val);
xlabel('D');
ylabel('H');
zlabel('M(H)');
shading interp;

figure(3);
surf(fd_max,fh_max,fsh_max);
shading interp;
grid on;
xlabel('D');
ylabel('H');
zlabel('\sigma(H)');



% coff_init2=[2.303289357,- 1.985958922717043,2.7264927140631485,0.18001531848872626];
% fun_sh=@(c1,c2,c3,c4,x,y) c1*y.^2.*x.^3+c2*y.^2.*x.^2+c3*y.*x+c4;
% sf2 = fit( [fd_val, fh_val], fsh_val, fun_sh,'StartPoint', coff_init2);

sf2 = fit( [fd_val, fh_val], fsh_val, 'poly23');

shv1=reshape(sf2([fd_val, fh_val]),61,length(H_list));
z=feval(sf2,[fd_val, fh_val]);
sum(abs(fsh_val-z));

params_s=coeffvalues(sf2);
figure(4);
plot(sf2,[fd_val, fh_val],fsh_val);
xlabel('D');
ylabel('H');
zlabel('\sigma(H)');
x=[params_m,params_s];

save(sprintf('%s/ML_init_params.mat',data_path),'fd_val','fh_val','fmh_val','fsh_val','-v7.3');
save(sprintf('%s/global_init_params.mat',data_path),'x','mhv1','shv1','fsh_val','-v7.3');