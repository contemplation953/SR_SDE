clc;clear;

root_path=pwd();
data_path=sprintf('%s/5-dof/data',root_path);
H_list=0.01:0.01:20;
d_list=0.005:0.0005:0.035;
c=load(sprintf('%s/global_params.mat',data_path)).x;

fd_max=zeros(61,length(H_list));
fh_max=zeros(61,length(H_list));
fmh_max=zeros(61,length(H_list));
fsh_max=zeros(61,length(H_list));

for i=1:61
    x=d_list(i);
    m_x=@(y) c(1)+c(2).*x+c(3).*y+c(4).*x.*y+c(5).*y.^2;
    s_x=@(y) c(6)+c(7).*x+c(8).*y+c(9).*x.^2+c(10).*x.*y+c(11).*y.^2+c(12).*x.^2.*y+c(13).*x.*y.^2+c(14).*y.^3;
    
    sim_m=m_x(H_list);
    sim_s=s_x(H_list);
    fmh_max(i,:)=sim_m';
    fsh_max(i,:)=sim_s';
    fd_max(i,:)=ones(1,length(H_list))*d_list(i);
    fh_max(i,:)=H_list;
end


%% 
figure(1);
surf(fd_max,fh_max,fmh_max);
grid on;
xlabel('D');
ylabel('H');
zlabel('M(H)');
shading interp;


figure(2);
surf(fd_max,fh_max,fsh_max);
grid on;
xlabel('D');
ylabel('H');
zlabel('\sigma^2(H)');
shading interp;
