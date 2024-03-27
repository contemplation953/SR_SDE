clc;clear;

root_path=pwd();
res=zeros(5,5);
oned=zeros(1,25);
index=1;
ssr_path=sprintf('%s/5-dof/ssr',root_path);
for i=1:5
    for j=1:5
        data=load(sprintf('%s/SSR_%d_%d.mat',ssr_path,i,j));
        res(i,j)=data.res;
        oned(index)=data.res;
        index=index+1;
    end
end


figure(1);
bar3(res,0.5,'grouped');
ylabel('d_m','FontSize',12);
zlabel('Loss_d','FontSize',12);
legend('d_{\sigma}=1','d_{\sigma}=2','d_{\sigma}=3','d_{\sigma}=4','d_{\sigma}=5','FontSize',12);
