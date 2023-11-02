clc;clear;

mh_p=4;
sh_p=5;
root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);
data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));

h=0.01;
H_list=0.01:h:25-h;
d_list=0.005:0.0005:0.035;
pdf_bound=1e-3;

fd_val=[];
fh_val=[];
fmh_val=[];
fsh_val=[];
for i=1:61
    data=load(sprintf('%s/ssr_%d_%d',ssr_path,mh_p,sh_p));
    coeff=data.x_list(i,:);
    coeff1=coeff(1:mh_p+1);
    coeff2=coeff(mh_p+2:end);
    [drift,diffusion]=FuncFactory(coeff1,coeff2);

    data2=load(sprintf('%s/H_%d.mat',data_path,i));
    monment_count=data2.monment_count;
    pdf=monment_count'/sum(monment_count)*200;
    idnex=find(pdf>=pdf_bound);
    begin_idnex=idnex(1);
    end_index=idnex(end);
    pst=H_list(begin_idnex:end_index);

    sim_m=drift(pst);
    sim_s=diffusion(pst);
    fmh_val=[fmh_val;sim_m'];
    fsh_val=[fsh_val;sim_s'];
    fd_val=[fd_val;ones(length(pst),1)*d_list(i)];
    fh_val=[fh_val;pst'];
end


%% 
figure(1);
plot3(fd_val,fh_val,fmh_val,'O-');
grid on;
xlabel('D');
ylabel('H');
zlabel('M(H)');
shading interp;
sf1=fit([fd_val, fh_val],fmh_val,'poly32');
params_m=coeffvalues(sf1);
figure(2);
plot(sf1,[fd_val, fh_val],fmh_val);
xlabel('D');
ylabel('H');
zlabel('M(H)');

figure(3);
plot3(fd_val,fh_val,fsh_val,'O-');
grid on;
xlabel('D');
ylabel('H');
zlabel('\sigma(H)');
shading interp;
sf2=fit([fd_val, fh_val],fsh_val,'poly32');
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
save(sprintf('%s/global_init_params.mat',data_path),'x','fsh_val','-v7.3');
