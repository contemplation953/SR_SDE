clc;clear;

mh_p=3;
sh_p=3;
root_path=pwd();
ssr_path=sprintf('%s/5-dof/ssr',root_path);
data_path=sprintf('%s/5-dof/data',root_path);
data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,mh_p,sh_p));

H_list=0.01:0.01:20;
d_list=0.005:0.0005:0.035;
pdf_bound=1e-3;

fd_val=[];
fh_val=[];
fmh_val=[];
fsh_val=[];

c=load(sprintf('%s/global_params.mat',data_path)).x;  
for i=1:61
    d=d_list(i);
    drift=@(h) -0.03071741433.*h - 0.011000170093282483.*d.*h.^2 + 0.02710054950 + 0.009704939708.*d.*h + 4.540426193.*d + 1.6259656446958923.*d.^2.*h + 5.584898066.*d.^2 + 2.*d.^3.*h;
    diffusion=@(h) 0.048377131.*h + 10.994370654887781.*d - d.^3.*h.^2 - 0.02269259174368238;

    data2=load(sprintf('%s/H_%d.mat',data_path,i));
    p_data=data2.p_data;
    idnex=find(p_data>=pdf_bound);
    begin_idnex=idnex(1);
    end_index=idnex(end);
    pst=H_list(begin_idnex:end_index);

    sim_m=drift(pst);
    sim_s=(diffusion(pst)).^2;
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


figure(2);
plot3(fd_val,fh_val,fsh_val,'O');
grid on;
xlabel('D');
ylabel('H');
zlabel('\sigma^2(H)');
shading interp;
