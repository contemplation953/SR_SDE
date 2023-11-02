clc;clear;

root_path=pwd();
res=zeros(5,5);
index=1;
oned=zeros(1,25);

for i=1:5
    for j=1:4
        data_file=sprintf('%s/N-DOF/10_dof_data/SSR_%d_%d.mat',root_path,i,j);
        if exist(data_file,'file')
            data=load(sprintf('%s/N-DOF/10_dof_data/SSR_%d_%d.mat',root_path,i,j));
            res(i,j)=data.res;
            oned(index)=data.res;
        else
            res(i,j)=0.2;
            oned(index)=0.2;
        end
        
        
        index=index+1;
    end
end

p=1:5;
x=repmat(p',1,5);
y=repmat(p,5,1);

figure(1);
bar3(res,0.5);
xlabel('The polynomial order of m(H)');
ylabel('The polynomial order of \sigma^2(H)');
zlabel('Error');

figure(2);
plot3(x,y,res,'-O');

figure(3);
plot(1:25,oned,'k*');